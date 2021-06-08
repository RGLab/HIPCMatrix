#' Process Affy
#'
#' @param input_files input file names
#' @param gef result of \code{ISCon$getDataset("gene_expression_files")} for one run.
#' @param metaData list of study-specific meta data
#'
#' @return background-corrected data.table of probe intensities
#'
#' @details
#' Requires annotation packages which can be installed from bioconductor.
#' \code{affy::justRMA()} will use the appropriate package.
#' Two are customCDF packages loaded from UpdateAnno: \code{huex10stv2cdf}
#' and \code{hursta2a520709cdf}
#'
.processAffy <- function(input_files, gef, metaData){
  # Background Correction Notes:
  # 'background' = TRUE performs function similar to normexp.fit.control and normexp.signal
  # from limma package.
  tmp <- getwd()
  setwd("/") # b/c filepaths are absolute and justRMA prepends wd
  eset <- affy::justRMA(filenames = input_files, normalize = FALSE, background = TRUE)
  setwd(tmp)
  exprs <- data.table(exprs(eset), keep.rownames = TRUE)
  setnames(exprs, "rn", "feature_id")

  # Names come from input_files. In case of isGeo, these are not exact
  # matches but usually have the gsm accession in them.
  if (any(grep("GSM", colnames(exprs)))) {
    nms <- grep("feature_id", colnames(exprs), invert = TRUE, value = TRUE)
    gsms <- regmatches(nms, regexpr("GSM\\d{6,7}", nms))
    setnames(exprs, nms, gsms)
  }

  return(exprs)
}
