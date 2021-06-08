#' Prep ImmPort Files
#'
#' @param study study accession eg \code{SDY269}
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
#' @param metaData list of study-specific meta data
#' @param input_files input file names
#'
#' @return path to raw, prepped input files
#'
.prepImmportFls <- function(study, gef, metaData, input_files){
  supp_files_dir <- .get_supp_files_dir(study, gef)
  if( metaData$platform == "Illumina") {
    mxList <- lapply(input_files, function(path){
      em <- fread(path)
      em <- .subsetIlluminaEM(em)
      em <- .prepIlluminaHeaders(em)
    })
  } else if (metaData$platform == "NA") {
    mxList <- lapply(input_files, fread)
    mxList <- .fixHeaders(mxList, study)
  }

  input_files <- .mxListToFlatFile(mxList, supp_files_dir, study)
  return(input_files)
}



#' Map experiment-sample or geo accessions to biosample accessions
#'
#' @param exprs data.table of expression
#'
.mapAccToBs <- function(exprs, gef){
  if (any(grepl("^(ES|GSM)\\d{6,7}$", colnames(exprs)))) {
    colToUse <- ifelse( any(grepl("ES", colnames(exprs))),
                        "expsample_accession",
                        "geo_accession")
  } else {
    colToUse <- "file_info_name"
  }
  nms <- grep("feature_id", colnames(exprs), value = TRUE, invert = TRUE)
  rep <- gef$biosample_accession[ match(nms, gef[[colToUse]]) ]

  # remove samples without matching biosample accession
  nms <- nms[!is.na(rep)]
  rep <- rep[!is.na(rep)]

  setnames(exprs, nms, rep)
  return(exprs)
}
