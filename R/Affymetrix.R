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
.process_affy <- function(input_files, gef, metaData){
  # Background Correction Notes:
  # 'background' = TRUE performs function similar to normexp.fit.control and normexp.signal
  # from limma package.
  wd <- getwd()
  setwd("/") # b/c filepaths are absolute and justRMA prepends wd
  eset <- affy::justRMA(filenames = input_files, normalize = FALSE, background = TRUE)
  setwd(wd)
  exprs_dt <- data.table(Biobase::exprs(eset), keep.rownames = TRUE)
  setnames(exprs_dt, "rn", "feature_id")

  # Names come from input_files. In case of isGeo, these are not exact
  # matches but usually have the gsm accession in them.
  if ( any(grep("GSM", colnames(exprs_dt))) ) {
    sample_names <- grep("feature_id",
                         colnames(exprs_dt),
                         invert = TRUE,
                         value = TRUE)
    gsms <- regmatches(sample_names, regexpr("GSM\\d{6,7}", sample_names))
    setnames(exprs_dt, sample_names, gsms)
  }

  return(exprs_dt)
}

