#' Normalize raw RNA-seq data
#'
#' @param counts_dt data.table with raw counts, one column per sample, one
#' row per feature
#'
#' @export
normalize_rnaseq <- function(counts_dt) {

  if (sum(is.na(exprs_dt)) > 0) {
    stop("Missing values found.")
  }

  feature_ids <- exprs_dt[ , feature_id ]
  exprs_dt[ , feature_id := NULL ]

  # Must ensure numeric values as conversion back to character can happen with
  # casting as matrix.
  exprs_mx <- as.matrix(apply(exprs_dt, 2, as.numeric))

  # newCountDataSet does not take duplicated column names, so assign temporary unique names
  orginal_colnames <- colnames(exprs_dt)
  colnames(exprs_dt) <- seq_len(ncol(exprs_dt))

  cds <- DESeq2::DESeqDataSetFromMatrix(countData = exprs_dt,
                                        colData = data.frame(cohort = rep(1, ncol(exprs_dt))),
                                        design = ~ 1)
  cds <- DESeq2::estimateSizeFactors(cds)
  cdsBlind <- DESeq2::estimateDispersions(cds)
  vsd <- DESeq2::varianceStabilizingTransformation(cdsBlind)
  norm_exprs <- SummarizedExperiment::assay(vsd)
  colnames(norm_exprs) <- orginal_colnames


  norm_exprs <- data.table(norm_exprs)
  norm_exprs[ , feature_id := rnames ]

  sample_id_columns <- grep("feature_id", names(norm_exprs), invert = TRUE)
  feature_id_column <- grep("feature_id", names(norm_exprs))
  setcolorder(norm_exprs, c(feature_id_column, sample_id_columns))

  norm_exprs

}


#' Normalize microarray data
#'
#' Normalize across samples to correct for technical variation using
#' quantile normalization. Also log2-transforms the data.
#'
#' @return matrix with normalized expression values.
#'
#' @param bg_corrected_dt data.table with background-corrected microarray
#' expression data
#' @param log2_transform (boolean) Should log2 transform be performed before normalization?
#' This should be FALSE for affy data which was read using RMA, which is already
#' in log2 scale.
#' @param force (boolean) Force log2 transform. If \code{log2_transform} is
#' \code{TRUE}, this function will error if \code{counts_dt} does not have
#' values greater than 100, unless \code{force = TRUE}.
#'
#' @export
normalize_microarray <- function(counts_dt,
                                 log2_transform = TRUE,
                                 force = FALSE) {

  if (sum(is.na(exprs_dt)) > 0) {
    stop("Missing values found.")
  }

  feature_ids <- exprs_dt[ , feature_id ]
  exprs_dt[ , feature_id := NULL ]

  # Must ensure numeric values as conversion back to character can happen with
  # casting as matrix.
  exprs_mx <- as.matrix(apply(exprs_dt, 2, as.numeric))
  cnames <- colnames(exprs_mx)

  # Do log2 transformation BEFORE normalization.
  if ( log2_transform ) {
    if ( max(exprs_mx) < 100 & !force) {
      stop("max(exprs_mx) < 100. ",
           "It is likely already in log2 scale. ",
           "Run with force=TRUE if you still want to log2 transform")
    }
    exprs_mx <- log2(exprs_mx)
  }

  norm_exprs <- preprocessCore::normalize.quantiles(exprs_mx)
  colnames(norm_exprs) <- cnames
  norm_exprs <- pmax(norm_exprs, 1)

  norm_exprs <- data.table(norm_exprs)
  norm_exprs[ , feature_id := rnames ]

  smpls <- grep("feature_id", names(norm_exprs), invert = TRUE)
  fid <- grep("feature_id", names(norm_exprs))
  setcolorder(norm_exprs, c(fid, smpls))

  norm_exprs
}


#' Normalize matrix
#'
#' @param exprs_dt data.table with raw, background-corrected expression. One
#' column per sample, one row per feature.
#' @param meta_data list of study-specific meta data
#'
normalize_matrix <- function(exprs_dt, platform) {
  switch(
    platform,
    "Illumina" = normalize_microarray(exprs_dt, log2_transform = TRUE),
    "Affymetrix" = normalize_microarray(exprs_dt, log2_transform = FALSE),
    "NA" = normalize_rnaseq(exprs_dt),
    "Stanford Functional Genomics Facility" = normalize_microarray(exprs_dt, log2_transform = FALSE),
    stop("Did not recognize platform: ", platform)
  )
}

