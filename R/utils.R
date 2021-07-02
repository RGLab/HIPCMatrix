#' Get supplemental files directory
#'
#' Create run-specific base directory for supplementary files.
#' Important that run specific directory is created to avoid overlap with previous runs.
#' This done with unique cell_type * arm_accession string (aka cohort_type)
#'
#' @param analysis_dir analysis directory
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
.get_supp_files_dir <- function(analysis_dir,
                                gef){
  supp_files_dir <- file.path(analysis_dir,
                              "supp_files",
                              paste0(unique(gef$type),
                                     "_",
                                     unique(gef$arm_accession)))
  dir.create(supp_files_dir, recursive = TRUE)
  return(supp_files_dir)
}

log_message <- function(msg) {
  message(sprintf("[%s] %s", Sys.time(), msg))
}

#' Write raw expression to a file
#'
#' Writes \code{ge_tbl} to a (tsv-formatted) \code{.txt} file, using controlled
#' format and file name.
#'
#' @return path to raw expression file
#'
#' @param ge_tbl object containing gene expression values (matrix, data.frame, etc)
#' @param supp_files_dir path to base directory (via \code{.get_supp_files_dir()})
#' @param study study accession eg \code{SDY269}
.write_raw_expression <- function(ge_tbl, supp_files_dir, study) {
  input_files <- file.path(supp_files_dir, paste0(study, "_raw_expression.txt"))
  fwrite(ge_tbl, file = input_files, sep = "\t", quote = FALSE, row.names = FALSE)
  return(input_files)
}

#' fix headers
#'
#' Custom fixes for individual studies
#'
#' @param ge_list list of matrices for one study
#' @param study study accession eg \code{SDY269}
.fix_headers <- function(ge_list, study){
  if ( study == "SDY224" ) {
    ge_list <- lapply(ge_list, function(x){
      setnames(x, as.character(x[1,]))
      x <- x[-(1:2),]
      colnames(x)[[1]] <- "ID_REF"
      return(x)
    })
  } else if ( study == "SDY400" ) {
    # Using mapping file provided by Hailong Meng at Yale, Dec 2018
    # since note in header of file is misleading due to gsm swaps made
    # later based on knowledge of switched samples.
    ge_list <- lapply(ge_list, function(x){
      mp <- fread("/share/files/Studies/SDY400/@files/rawdata/gene_expression/SDY400_HeaderMapping.csv")
      setnames(x, colnames(x), as.character(x[2,]))
      x <- x[-(1:2),]
      smpls <- grep("SAMPLE", colnames(x), value = T)
      titles <- mp$Title[ match(smpls, mp$Sample) ]
      setnames(x, smpls, titles)
      return(x)
    })
  } else if ( study == "SDY1325" ) {
    ge_list <- lapply(ge_list, function(x){
      setnames(x, colnames(x), as.character(x[5,]))
      x <- x[6:nrow(x),]
      return(x)
    })
  } else if ( study == "SDY1324" ) {
    # Custom header mapping provided by authors via P.Dunn Dec 2018.
    ge_list <- lapply(ge_list, function(x){
      mp <- fread("/share/files/Studies/SDY1324/@files/rawdata/gene_expression/raw_counts/SDY1324_Header_Mapping.csv")
      accs <- grep("V1", colnames(x), invert = TRUE, value = TRUE)
      esNms <- mp$experimentAccession[ match(accs, mp$AuthorGivenId) ]
      setnames(x, accs[!is.na(esNms)], esNms[!is.na(esNms)])
      return(x)
    })
  } else if (study == "SDY787") {
    # Fist number is unique id
    ge_list <- lapply(ge_list, function(x) {
      # Remove first "_" and everything following
      setnames(x, colnames(x), gsub("_.*$", "", colnames(x)))
    })
  }
  return(ge_list)
}

#' matrix list to flat file
#'
#' ge_list to flat file. Used when samples are in different files on GEO
#'
#' @return path to raw, prepped input files
#'
#' @param ge_list list of gene expression tables (matrix, data.frame, etc)
#' @param supp_files_dir path to base directory (via \code{.get_supp_files_dir()})
#' @param study study accession eg \code{SDY269}
.ge_list_to_flat_file <- function(ge_list, supp_files_dir, study){
  ge_df <- Reduce(f = function(x, y) {merge(x, y)}, ge_list)
  input_files <- .write_raw_expression(ge_df, supp_files_dir, study)
}
#######################################
###            MAPPING              ###
#######################################



#' Summarize Matrix
#'
#' Summarize any duplicated genes so that there is
#' one entry per gene.
#'
#' @details This function summarizes normalized expression values by mean
#' into a table with one entry per gene name based on \code{feature_gene_map}.
#' The input data should be in log2 space.
#'
#' @param ge_dt data.table of gene expression values. Features should be listed in
#' \code{feature_id} column, with additional columns for each sample.
#' @param feature_gene_map data.table mapping original gene id (eg probe id, ensembl id, gene alias)
#' to gene symbol. It should have two columns: \code{featureid} (original) and \code{genesymbol} (updated).
#'
#' @export
summarize_by_gene_symbol <- function(ge_dt,
                             feature_gene_map){
  em <- copy(ge_dt)
  em[ , gene_symbol := feature_gene_map[match(em$feature_id, feature_gene_map$featureid), genesymbol] ]
  em <- em[ !is.na(gene_symbol) & gene_symbol != "NA" ]
  sum_exprs <- em[ , lapply(.SD, mean),
                   by = "gene_symbol",
                   .SDcols = grep("^BS", colnames(em)) ]
}

#' Write matrix to file system
#'
#' @details Writes four versions of flat tsv files:
#'   1. \code{<matrix_name>.raw.tsv}: raw, background-corrected values
#'   2. \code{<matrix_name>.tsv}: normalized values
#'   3. \code{<matrix_name>.summary.tsv}: normalized values, summarized by gene symbol
#'   (based on current annotation)
#'   4. \code{<matrix_name>.summary.orig}: normalized values, summarized by gene symbol
#'   (based on original annotation)
#'
#' @param pipeline.root pipeline.root
#' @param matrix_name Name of matrix
#' @param exprs data.table of raw (background-corrected) expression (from \code{makeRawMatrix()})
#' @param norm_exprs data.table of normalized expression (from \code{normalizeMatrix()})
#' @param sum_exprs data.table of summarized expression (from \code{summarizeMatrix()})
#' @param verbose Write verbose print statements?
#'
#' @export
write_matrix <- function(output_dir,
                         matrix_name,
                         exprs,
                         norm_exprs,
                         sum_exprs,
                         verbose = FALSE){

  .write_em <- function(df, file_name, verbose){
    if ( verbose ) message("Writing ", file_name, "...")
    fwrite(df,
           file = file.path(output_dir,
                            file_name),
           sep = "\t",
           quote = FALSE,
           row.names = FALSE)
  }

  # Raw EM
  .write_em(exprs, paste0(matrix_name, ".tsv.raw"), verbose)

  # Normalized EM
  .write_em(norm_exprs, paste0(matrix_name, ".tsv"), verbose)

  # summary EM
  .write_em(sum_exprs, paste0(matrix_name, ".tsv.summary"), verbose)

  # original summary EM assuming run created with _orig FasId
  .write_em(sum_exprs, paste0(matrix_name, ".tsv.summary.orig"), verbose)

}


parse_pipeline_inputs <- function(pipeline_inputs) {
  list(
    study = gsub("/Studies/", "", pipeline_inputs$labkey.url.path),
    matrix_name = gsub("\\.tsv", "", pipeline_inputs$output.tsv),
    selected_biosamples = pipeline_inputs$selectedBiosamples,
    fas_id = pipeline_inputs$fasId,
    labkey.url.base = pipeline_inputs$labkey.url.base,
    base_dir = pipeline_inputs$pipeline.root,
    taskOutputParams = pipeline_inputs$taskOutputParams,
    verbose = TRUE
  )
}


get_input_params <- function(con = CreateConnection(""),
                             matrix_name) {
  # First make sure the matrix exists
  stopifnot(matrix_name %in% con$listGEMatrices()$name)

  # First check for inputParams
  baseUrl <- con$config$labkey.url.base
  study <- con$listGEMatrices()[name == matrix_name, folder]
  if ( Rlabkey::labkey.webdav.pathExists(
    baseUrl = con$config$labkey.url.base,
    folderPath = paste0("Studies/", study),
    remoteFilePath = file.path("rawdata",
                               "gene_expression",
                               "create-matrix",
                               matrix_name,
                               "create-matrix-vars.tsv"),
    fileSet='@files'
  ) ) {
    local_path <- tempfile()
    Rlabkey::labkey.webdav.get(baseUrl = con$config$labkey.url.base,
                               folderPath = paste0("Studies/", study),
                               remoteFilePath = file.path("rawdata",
                                                          "gene_expression",
                                                          "create-matrix",
                                                          matrix_name,
                                                          "create-matrix-vars.tsv"),
                               localFilePath = local_path,
                               fileSet='@files')
    input_params <- fread(local_path)
    return(parse_pipeline_inputs(input_params))
  } else {
    # Get biosamples and fasID
    input_samples <- Rlabkey::labkey.selectRows(
      baseUrl = con$config$labkey.url.base,
      folderPath = paste0("Studies/", study),
      schemaName = "assay.ExpressionMatrix.matrix",
      queryName = "inputSamples_computed",
      colSelect = c("Biosample", "Run/featureSet"),
      colFilter = makeFilter(c("Run/Name", "EQUAL", matrix_name)),
      colNameOpt = "fieldname"
    )
    selected_biosamples <- paste0(input_samples$Biosample, collapse = ",")
    fas_id <- unique(input_samples$`Run/featureSet`)
    return(list(
      study = study,
      matrix_name = matrix_name,
      selected_biosamples = selected_biosamples,
      fas_id = fas_id,
      labkey.url.base = baseUrl,
      base_dir = file.path("/share", "files", "Studies", study, "@files"),
      verbose = TRUE
    ))
  }
}
