## Helpers

#' Function to update annotation of a character vector of gene aliases
#'
#' @param aliases character vector of gene aliases
#' @import data.table
#' @export
#'
mapAlias2Symbol <- function(aliases){
  vec_dt <- data.table(aliases)
  setnames(vec_dt, "aliases", "ALIAS")
  vec_dt[hgncAlias2Symbol, SYMBOL := SYMBOL, on = c(ALIAS = "ALIAS")]
  return(vec_dt$SYMBOL)
}


#' Function to map probe ids to gene symbols from an annotation package
#'
#' @param annoPkg name of annotation package to use
#' @param outPath file path for tsv to write out for upload to ImmuneSpace
#' @export
#'
makeAnnoDf <- function(annoPkg, outPath = NULL){
  library(annoPkg, character.only = TRUE)
  prbLs <- paste0(gsub("\\.db", "", annoPkg), "ENTREZID")
  probes <- ls(eval(parse(text = prbLs)))
  symLs <- gsub("ENTREZID", "SYMBOL", prbLs)
  syms <- unlist(mget(probes, eval(parse(text = symLs))))
  syms <- syms[ !is.na(syms)] # No NAs or "" allowed in IS FAS
  res <- data.frame(Probe_ID = names(syms),
                    Gene_Symbol = syms,
                    stringsAsFactors = F)
  if( !is.null(outPath) ){
    write.table(res,
                file = outPath,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE)
  }
  return(res)
}



#' @export updateFAS
updateFAS <- function(baseUrl, folderPath = "/Studies/", fasNms = NULL){

  # vars ---------------------------------------------------
  schemaName <- "Microarray"

  # helper fn ----------------------------------------------
  getAnno <- function(nm, currFAS, baseUrl){
    annoSetId <- currFAS$RowId[ currFAS$Name == nm ]
    fasQuery <- sprintf("SELECT * from FeatureAnnotation
                          where FeatureAnnotationSetId='%s';",
                        annoSetId)
    features <- labkey.executeSql(baseUrl = baseUrl,
                                  folderPath = folderPath,
                                  schemaName = schemaName,
                                  sql = fasQuery,
                                  colNameOpt = "fieldname",
                                  showHidden = TRUE)
    toDrop <- c("Created", "CreatedBy", "Modified", "ModifiedBy")
    features <- features[ , !(colnames(features) %in% toDrop) ]
  }

  currFas <- function(baseUrl){
    res <- data.table(labkey.selectRows(baseUrl = baseUrl,
                                        folderPath = folderPath,
                                        schemaName = schemaName,
                                        queryName = "FeatureAnnotationSet",
                                        colNameOpt = "fieldname",
                                        showHidden = TRUE ))
  }


  # MAIN ------------------------------------------------------------
  # for each name in fas$name see if there is an updated version then work on the
  # updated version if there is one or create one if there is not. Exceptions are those
  # with "Do Not Update", "non-updatable" or similar string in the Comment column.
  currFAS <- currFas(baseUrl)

  if (is.null(fasNms)) {
    fasNms <- currFAS$Name[ grep("[N|n]o[n|t|][-| ][u|U]pdat[e|able]",
                                 currFAS$Comment, invert = TRUE) ]
  }

  lapply(fasNms, FUN = function(nm){
    message(paste0("Updating: ", nm))
    currAnno <- getAnno(nm, currFAS, baseUrl)
    orNm <- paste0(nm, "_orig")

    if( orNm %in% currFAS$Name ){
      message("Updating FAS")
      # if orig is present means that update has been performed at least once.
      # Update the rows of the previously updated anno using the orig
      # as the base for mapping to ensure that updates of updates are avoided
      # and the correct rowIds remain.
      orAnno <- getAnno(orNm, currFAS, baseUrl)
      orAnno$GeneSymbol <- gsub("---", "", orAnno$GeneSymbol)
      orAnno$GeneSymbol <- mapAlias2Symbol(orAnno$GeneSymbol)
      currAnno$GeneSymbol <- orAnno$GeneSymbol[ match(currAnno$FeatureId, orAnno$FeatureId) ]
      currAnno[ is.na(currAnno) ] <- ""
      toUpdate <- data.frame(currAnno, stringsAsFactors = FALSE)
      done <- labkey.updateRows(baseUrl = baseUrl,
                                folderPath = folderPath,
                                schemaName = schemaName,
                                queryName = "FeatureAnnotation",
                                toUpdate = toUpdate)
    }else{
      print("Creating New updated FAS")
      # First create featureAnnotationSet with "_orig" name and original annotation
      toImport <- data.frame(currFAS[ currFAS$Name == nm, ])
      toImport <- toImport[ , !(colnames(toImport) == "RowId") ]
      toImport$Name <- orNm
      toImport$Comment <- "Do not update"
      addFas <- labkey.importRows(baseUrl = baseUrl,
                                  folderPath = folderPath,
                                  schemaName = schemaName,
                                  queryName = "FeatureAnnotationSet",
                                  toImport = toImport)

      # check that new "_orig" set has been imported correctly
      nowFas <- currFas(baseUrl)
      if( !(toImport$Name[[1]] %in% nowFas$Name) ){
        stop("Original FAS (",toImport$Name[[1]],") not imported correctly")
      }

      # Prep for importRows by removing rowIds (will be given) and using new FASid
      orAnno <- currAnno[ , !(colnames(currAnno) == "RowId") ]
      orAnno$FeatureAnnotationSetId <- nowFas$RowId[ nowFas$Name == toImport$Name[[1]] ]
      orAnno[ is.na(orAnno) ] <- ""
      toImport <- data.frame(orAnno, stringsAsFactors = FALSE)
      addFeatures <- labkey.importRows(baseUrl = baseUrl,
                                       folderPath = folderPath,
                                       schemaName = schemaName,
                                       queryName = "FeatureAnnotation",
                                       toImport = toImport)

      featureChk <- labkey.selectRows(baseUrl = baseUrl,
                                      folderPath = folderPath,
                                      schemaName = schemaName,
                                      queryName = "FeatureAnnotation",
                                      colNameOpt = "fieldname",
                                      colSelect = c("GeneSymbol", "FeatureId"),
                                      colFilter = makeFilter(c("FeatureAnnotationSetId",
                                                               "EQUALS",
                                                               unique(toImport$FeatureAnnotationSetId)))
      )
      featureChk[ is.na(featureChk) ] <- ""
      featureChk <- featureChk[ order(featureChk$FeatureId), ]
      imported <- data.frame(GeneSymbol = toImport$GeneSymbol,
                             FeatureId = toImport$FeatureId,
                             stringsAsFactors = FALSE)
      imported <- imported[ order(imported$FeatureId), ]
      rownames(imported) <- rownames(featureChk) <- NULL

      if( all.equal(featureChk, imported) ){
        # Now update the old fasId rows with new geneSymbols
        currAnno$GeneSymbol <- mapAlias2Symbol(currAnno$GeneSymbol)
        currAnno[ is.na(currAnno) ] <- ""
        FAUpdate <- data.frame(currAnno, stringsAsFactors = FALSE)
        FAdone <- labkey.updateRows(baseUrl = baseUrl,
                                    folderPath = folderPath,
                                    schemaName = schemaName,
                                    queryName = "FeatureAnnotation",
                                    toUpdate = FAUpdate)

      }else{
        stop("Original FA not uploaded correctly to *_orig table")
      }
    }

    # Update FAS$comment to be packageVersion of org.Hs.eg.db
    updateFAS <- data.frame(currFAS[ currFAS$Name == nm, ], stringsAsFactors = FALSE)
    updateFAS$Comment <- paste0("Alias2Symbol mapping with HGNC dataset from: ",
                                UpdateAnno::hgncAlias2Symbol_version)
    FASdone <- labkey.updateRows(baseUrl = baseUrl,
                                 folderPath = folderPath,
                                 schemaName = schemaName,
                                 queryName = "FeatureAnnotationSet",
                                 toUpdate = updateFAS)
  })
  return(TRUE)
}

#######################################################################
###            2. Update the expression matrices                    ###
#######################################################################
# if running by hand on rsT/rsP
# library(ImmuneSpaceR)
# library(Rlabkey)
# library(data.table)


#' @export updateEMs
updateEMs <- function(sdy, runsDF, baseUrl, folderPath = "/Studies/"){
  print(paste0("working on study: ", sdy))

  # get file basenames present on server
  sdyPath <- ifelse( folderPath == "/Studies/",
                     paste0(folderPath, sdy),
                     folderPath)
  dirPath <- paste0("/share/files", sdyPath, "/@files/analysis/exprs_matrices")
  fls <- list.files(dirPath)
  tmp <- unique(unlist(strsplit(fls, split = ".tsv", fixed = TRUE)))
  baseNms <- tmp[ !(tmp %in% c(".summary",".summary.orig", ".raw", ".immsig")) ]

  # checks
  if(!all(baseNms %in% runsDF$name)){
    baseNms <- baseNms[ baseNms %in% runsDF$name]
    warning("Extra files / basenames present in current study.  Please delete.")
  }

  if( !(all(runsDF$name[ runsDF$folder_name == sdy] %in% baseNms))){
    stop("Runs missing from files!")
  }

  # go through each baseNm to update summary tsv
  sapply(baseNms, function(nm){
    print(paste0("working on baseNm: ", nm))
    baseFls <- fls[ grep(nm, fls) ]

    # Rename original summary file to tsv.summary.orig if necessary (first time only)
    # As of 7/2018, runCreateMatrix() is creating .summary.orig file at time of first run
    if( !(paste0(nm, ".tsv.summary.orig") %in% baseFls) ){
      sumFl <- paste0(dirPath, "/", nm, ".tsv.summary")
      dmp <- file.rename( sumFl, paste0(sumFl, ".orig") )
    }

    # get probe-level original df and update annotation using only features for FASid
    # Reason for limiting to FASid and therefore needing executeSql instead of SelectRows
    # is that original FAS entries may have same probes mapped to different genes based
    # on changes in bioconductor libraries over time.  ExecuteSql is only way to get FASid
    # because it is a lookup even though it is in microarray.FeatureAnnotation.
    prbEM <- fread(file.path(dirPath, paste0(nm, ".tsv")))
    annoSetId <- runsDF$featureset[ runsDF$name == nm ]
    currFAS <- data.table(labkey.selectRows(baseUrl = baseUrl,
                                            folderPath = "/Studies/",
                                            schemaName = "microarray",
                                            queryName = "FeatureAnnotationSet",
                                            colNameOpt = "fieldname",
                                            showHidden = TRUE ))
    currAnnoNm <- currFAS$Name[ currFAS$RowId == annoSetId]

    # Handle possibility that _orig anno was used to create mx (e.g. annotation
    # updated before matrix created.)
    if(grepl("_orig", currAnnoNm)){
      annoSetId <- currFAS$RowId[ currFAS$Name == gsub("_orig","", currAnnoNm)]
    }
    sqlStr <- sprintf("SELECT FeatureAnnotationSetId, FeatureId, GeneSymbol
                    from FeatureAnnotation
                    where FeatureAnnotationSetId='%s';", annoSetId)
    features <- labkey.executeSql(baseUrl = baseUrl,
                                  folderPath = "/Studies/",
                                  schemaName = "Microarray",
                                  sql = sqlStr,
                                  colNameOpt = "fieldname")
    prbEM[ , feature_id := as.character(feature_id)] # for SDY80 where probes have integer vals
    prbEM <- prbEM[features, gene_symbol := GeneSymbol, on = c(feature_id = "FeatureId")]

    # Summarize - lifted from Create-Matrix.R
    em <- prbEM[ !is.na(gene_symbol) & gene_symbol != "NA" ]
    sumEM <- em[ , lapply(.SD, mean), by="gene_symbol", .SDcols = grep("^BS", colnames(em)) ]

    write.table(sumEM, file = paste0(dirPath, "/", nm, ".tsv.summary"),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
  })
  return(TRUE)
}
