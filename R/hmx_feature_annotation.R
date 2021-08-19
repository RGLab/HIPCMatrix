#' @include HMX.R
NULL

## Helpers

#' map alias to gene symbol
#'
#' Function to update annotation of a character vector of gene aliases, using
#' the mapping provided by \link{\code{hgncAlias2Symbol}}
#'
#' @param aliases character vector of gene aliases
#' @import data.table
#' @export
#'
mapAlias2Symbol <- function(aliases,
                            gene_symbol_map = hgncAlias2Symbol) {
  vec_dt <- data.table(aliases)
  setnames(vec_dt, "aliases", "ALIAS")
  vec_dt[gene_symbol_map, SYMBOL := SYMBOL, on = c(ALIAS = "ALIAS")]
  return(vec_dt$SYMBOL)
}


#' make annotation df
#'
#' Function to map probe ids to gene symbols from an annotation package
#'
#' @param annoPkg name of annotation package to use
#' @param outPath file path for tsv to write out for upload to ImmuneSpace
#' @export
#'
makeAnnoDf <- function(annoPkg, outPath = NULL) {
  library(annoPkg, character.only = TRUE)
  prbLs <- paste0(gsub("\\.db", "", annoPkg), "ENTREZID")
  probes <- ls(eval(parse(text = prbLs)))
  symLs <- gsub("ENTREZID", "SYMBOL", prbLs)
  syms <- unlist(mget(probes, eval(parse(text = symLs))))
  syms <- syms[!is.na(syms)] # No NAs or "" allowed in IS FAS
  res <- data.frame(
    Probe_ID = names(syms),
    Gene_Symbol = syms,
    stringsAsFactors = F
  )
  if (!is.null(outPath)) {
    write.table(res,
      file = outPath,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )
  }
  return(res)
}

##### ----------------- HMX methods ------------------------ #####

# Update feature annotation set tables.
# see HMX$updateFAS.
updateFAS <- function(con, fas_names = NULL) {
  if (grepl("HIPC", con$config$labkey.url.path)) {
    stop("Don't update ImmuneSignatures annotation!")
  }

  # vars ---------------------------------------------------
  schemaName <- "Microarray"

  # helper fn ----------------------------------------------
  getAnno <- function(fas_name, currFAS) {
    annoSetId <- currFAS$RowId[currFAS$Name == fas_name]
    fasQuery <- sprintf(
      "SELECT * from FeatureAnnotation
                          where FeatureAnnotationSetId='%s';",
      annoSetId
    )
    features <- labkey.executeSql(
      baseUrl = con$config$labkey.url.base,
      folderPath = con$config$labkey.url.path,
      schemaName = schemaName,
      sql = fasQuery,
      colNameOpt = "fieldname",
      showHidden = TRUE
    )
    toDrop <- c("Created", "CreatedBy", "Modified", "ModifiedBy")
    features <- features[, !(colnames(features) %in% toDrop)]
  }

  # MAIN ------------------------------------------------------------
  # for each name in fas$name see if there is an updated version then work on the
  # updated version if there is one or create one if there is not. Exceptions are those
  # with "Do Not Update", "non-updatable" or similar string in the Comment column.
  currFAS <- ImmuneSpaceR:::.getLKtbl(
    con,
    schema = schemaName,
    query = "FeatureAnnotationSet",
    colNameOpt = "fieldname"
  )

  if (is.null(fas_names)) {
    # If single study, only update FAS associated with matrices in that study
    if (con$study == "Studies") {
      fas_names <- currFAS$Name[grep("[N|n]o[n|t|][-| ][u|U]pdat[e|able]",
        currFAS$Comment,
        invert = TRUE
      )]
    } else {
      sdy_fas <- currFAS[RowId %in% unique(con$listGEMatrices()$featureset)]$Name
      sdy_fas <- gsub("_orig", "", sdy_fas)
      fas_names <- currFAS[!grepl(
        "[N|n]o[n|t|][-| ][u|U]pdat[e|able]",
        currFAS$Comment
      ) &
        currFAS$Name == sdy_fas]$Name
      if (length(fas_names) < 1) {
        stop("Could not determine fas_name for this study. Please specify fas_names")
      }
    }
  }

  lapply(fas_names, FUN = function(fas_name) {
    log_message(paste0("Updating: ", fas_name))
    curr_anno <- getAnno(fas_name, currFAS)
    orig_name <- paste0(fas_name, "_orig")

    if (orig_name %in% currFAS$Name) {
      log_message("Updating FAS")
      # if orig is present means that update has been performed at least once.
      # Update the rows of the previously updated anno using the orig
      # as the base for mapping to ensure that updates of updates are avoided
      # and the correct rowIds remain.
      orig_anno <- getAnno(orig_name, currFAS)
      orig_anno$GeneSymbol <- gsub("---", "", orig_anno$GeneSymbol)
      orig_anno$GeneSymbol <- mapAlias2Symbol(orig_anno$GeneSymbol)
      curr_anno$GeneSymbol <- orig_anno$GeneSymbol[match(curr_anno$FeatureId, orig_anno$FeatureId)]
      curr_anno[is.na(curr_anno)] <- ""
      toUpdate <- data.frame(curr_anno, stringsAsFactors = FALSE)
      done <- labkey.updateRows(
        baseUrl = con$config$labkey.url.base,
        folderPath = "/Studies/",
        schemaName = schemaName,
        queryName = "FeatureAnnotation",
        toUpdate = toUpdate
      )
    } else {
      print("Creating New updated FAS")
      # First create featureAnnotationSet with "_orig" name and original annotation
      toImport <- data.frame(currFAS[currFAS$Name == fas_name, ])
      toImport <- toImport[, !(colnames(toImport) == "RowId")]
      toImport$Name <- orig_name
      toImport$Comment <- "Do not update"
      addFas <- labkey.importRows(
        baseUrl = con$config$labkey.url.base,
        folderPath = "/Studies/",
        schemaName = schemaName,
        queryName = "FeatureAnnotationSet",
        toImport = toImport
      )

      # check that new "_orig" set has been imported correctly
      nowFas <- ImmuneSpaceR:::.getLKtbl(
        con,
        schema = schemaName,
        query = "FeatureAnnotationSet",
        colNameOpt = "fieldname"
      )
      if (!(toImport$Name[[1]] %in% nowFas$Name)) {
        stop("Original FAS (", toImport$Name[[1]], ") not imported correctly")
      }

      # Prep for importRows by removing rowIds (will be given) and using new FASid
      orig_anno <- curr_anno[, !(colnames(curr_anno) == "RowId")]
      orig_anno$FeatureAnnotationSetId <- nowFas$RowId[nowFas$Name == toImport$Name[[1]]]
      orig_anno[is.na(orig_anno)] <- ""
      toImport <- data.frame(orig_anno, stringsAsFactors = FALSE)
      addFeatures <- labkey.importRows(
        baseUrl = con$config$labkey.url.base,
        folderPath = "/Studies/",
        schemaName = schemaName,
        queryName = "FeatureAnnotation",
        toImport = toImport
      )

      featureChk <- labkey.selectRows(
        baseUrl = con$config$labkey.url.base,
        folderPath = "/Studies/",
        schemaName = schemaName,
        queryName = "FeatureAnnotation",
        colNameOpt = "fieldname",
        colSelect = c("GeneSymbol", "FeatureId"),
        colFilter = makeFilter(c(
          "FeatureAnnotationSetId",
          "EQUALS",
          unique(toImport$FeatureAnnotationSetId)
        ))
      )
      featureChk[is.na(featureChk)] <- ""
      featureChk <- featureChk[order(featureChk$FeatureId), ]
      imported <- data.frame(
        GeneSymbol = toImport$GeneSymbol,
        FeatureId = toImport$FeatureId,
        stringsAsFactors = FALSE
      )
      imported <- imported[order(imported$FeatureId), ]
      rownames(imported) <- rownames(featureChk) <- NULL

      if (all.equal(featureChk, imported)) {
        # Now update the old fasId rows with new geneSymbols
        curr_anno$GeneSymbol <- mapAlias2Symbol(curr_anno$GeneSymbol)
        curr_anno[is.na(curr_anno)] <- ""
        FAUpdate <- data.frame(curr_anno, stringsAsFactors = FALSE)
        FAdone <- labkey.updateRows(
          baseUrl = con$config$labkey.url.base,
          folderPath = "/Studies/",
          schemaName = schemaName,
          queryName = "FeatureAnnotation",
          toUpdate = FAUpdate
        )
      } else {
        stop("Original FA not uploaded correctly to *_orig table")
      }
    }

    # Update FAS$comment to indicate date of gene map update from HUGO
    updateFAS <- data.frame(currFAS[currFAS$Name == fas_name, ], stringsAsFactors = FALSE)
    updateFAS$Comment <- paste0(
      "Alias2Symbol mapping with HGNC dataset from: ",
      hgncAlias2Symbol_version
    )
    FASdone <- labkey.updateRows(
      baseUrl = con$config$labkey.url.base,
      folderPath = "/Studies/",
      schemaName = schemaName,
      queryName = "FeatureAnnotationSet",
      toUpdate = updateFAS
    )
  })
  return(TRUE)
}

# Update expression matrices
# See HMX$updateEMs
updateEMs <- function(con) {
  if (!grepl("SDY", con$study)) {
    stop("updateEMs must be run one study at a time.")
  }

  log_message("Updating annotation for ", con$study)

  ge_dir <- file.path(
    "/share/files/Studies", con$study,
    "@files/analysis/exprs_matrices"
  )
  ge_files <- list.files(ge_dir)
  tmp <- unique(unlist(strsplit(ge_files, split = ".tsv", fixed = TRUE)))
  base_names <- tmp[!(tmp %in% c(".summary", ".summary.orig", ".raw", ".immsig"))]

  if (!all(base_names %in% con$listGEMatrices())) {
    base_names <- base_names[base_names %in% con$listGEMatrices()]
    warning("Extra files / basenames present in current study.  Please delete.")
  }

  if (!(all(con$listGEMatrices()$name %in% base_names))) {
    stop("Runs missing from files!")
  }

  # Go through each base_name to update summary tsv
  sapply(base_names, function(mx_name) {
    log_message("Updating ", mx_name)
    mx_files <- ge_files[grep(mx_name, ge_files)]


    # Rename original summary file to tsv.summary.orig if necessary (first time only)
    # For legacy use only, as HIPCMatrix pipeline should create .summary.orig
    # file with first run.
    if (!(paste0(mx_name, ".tsv.summary.orig") %in% mx_files)) {
      sumFl <- paste0(ge_dir, "/", mx_name, ".tsv.summary")
      dmp <- file.rename(sumFl, paste0(sumFl, ".orig"))
    }

    # get probe-level original df and update annotation using only features for FASid
    # Reason for limiting to FASid and therefore needing executeSql instead of SelectRows
    # is that original FAS entries may have same probes mapped to different genes based
    # on changes in bioconductor libraries over time.  ExecuteSql is only way to get FASid
    # because it is a lookup even though it is in microarray.FeatureAnnotation.
    prb_mx <- fread(file.path(ge_dir, paste0(mx_name, ".tsv")))
    anno_set_id <- con$listGEMatrices()$featureset[con$listGEMatrices()$name == mx_name]
    curr_fas <- data.table(labkey.selectRows(
      baseUrl = con$config$labkey.url.base,
      folderPath = "/Studies/",
      schemaName = "microarray",
      queryName = "FeatureAnnotationSet",
      colNameOpt = "fieldname",
      showHidden = TRUE
    ))
    curr_anno_name <- curr_fas$Name[curr_fas$RowId == anno_set_id]


    # Handle possibility that _orig anno was used to create mx (e.g. annotation
    # updated before matrix created.)
    if (grepl("_orig", curr_anno_name)) {
      anno_set_id <- curr_fas$RowId[curr_fas$Name == gsub("_orig", "", curr_anno_name)]
    }


    sqlStr <- sprintf("SELECT FeatureAnnotationSetId, FeatureId, GeneSymbol
                    from FeatureAnnotation
                    where FeatureAnnotationSetId='%s';", anno_set_id)
    features <- data.table(labkey.executeSql(
      baseUrl = con$config$labkey.url.base,
      folderPath = "/Studies/",
      schemaName = "Microarray",
      sql = sqlStr,
      colNameOpt = "rname"
    ))

    prb_mx[, feature_id := as.character(feature_id)] # for SDY80 where probes have integer vals
    sum_mx <- summarize_by_gene_symbol(
      prb_mx,
      feature_gene_map = features,
      method = "mean",
      verbose = con$config$verbose
    )

    fwrite(
      sum_mx,
      file = file.path(
        ge_dir,
        paste0(mx_name, ".tsv.summary")
      ),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  })
  invisible(con)
}
