#' Normalize raw RNA-seq data
#'
#' RNA data will be normalized using the variance stabilizing transformation,
#' as implemented by the \code{DESeq2::vst} function.
#'
#' @return a matrix with one row per feature and one column per sample.
#'
#' @param counts_mx matrix with raw counts, one column per sample, one
#' row per feature
#'
#' @export
normalize_rnaseq <- function(counts_mx,
                             verbose = FALSE) {

  if (verbose) message(" --- normalize_rnaseq --- ")
  if (verbose) message("Normalizing counts data using variance stabilizing transformation...")
  if (sum(is.na(counts_mx)) > 0) {
    stop("Missing values found.")
  }
  if ( sum(duplicated(colnames(counts_mx))) > 0 ) {
    warning("Duplicate column name: ", colnames(counts_mx)[duplicated(colnames(counts_mx))])
  }
  # newCountDataSet does not take duplicated column names, so assign temporary unique names
  original_colnames <- colnames(counts_mx)
  colnames(counts_mx) <- seq_len(ncol(counts_mx))

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_mx,
                                        colData = data.frame(sample = original_colnames),
                                        design = ~ 1)
  vsd <- DESeq2::vst(dds)

  norm_exprs <- SummarizedExperiment::assay(vsd)
  colnames(norm_exprs) <- original_colnames

  norm_exprs
}


#' Normalize microarray data
#'
#' Normalize across samples to correct for technical variation using
#' quantile normalization. Also log2-transforms the data.
#'
#' @return a matrix with one row per feature and one column per sample.
#'
#'
#' @param exprs_mx matrix with background-corrected microarray
#' expression data, one column per sample, one row per feature.
#' @param log2_transform (boolean) Should log2 transform be performed before normalization?
#' This should be FALSE for affy data which was read using RMA, which is already
#' in log2 scale.
#' @param force (boolean) Force log2 transform. If \code{log2_transform} is
#' \code{TRUE}, this function will error if \code{counts_dt} does not have
#' values greater than 100, unless \code{force = TRUE}.
#'
#' @export
normalize_microarray <- function(exprs_mx,
                                 log2_transform = TRUE,
                                 force = FALSE,
                                 verbose = FALSE) {

  if (verbose) message(" --- normalize_microarray --- ")
  if (sum(is.na(exprs_mx)) > 0) {
    stop("Missing values found.")
  }
  if ( sum(duplicated(colnames(exprs_mx))) > 0 ) {
    warning("Duplicate column name: ", colnames(exprs_mx)[duplicated(colnames(exprs_mx))])
  }
  # normalize.quantiles removes row and column names
  cnames <- colnames(exprs_mx)
  rnames <- rownames(exprs_mx)

  # Do log2 transformation BEFORE normalization.
  if ( log2_transform ) {
    if ( max(exprs_mx) < 100 )
      if ( !force ) {
        stop("max(exprs_mx) < 100. ",
             "It is likely already in log2 scale. ",
             "Run with force=TRUE if you still want to log2 transform")
      } else if ( verbose ) {
        message("max(exprs_mx) < 100. Forcing log2 transform... ")
      }
    if (verbose) message("log2-transforming exprs_mx")
    exprs_mx <- log2(exprs_mx + 1)
  }
  if (verbose) message("Performing quantile normalization...")
  norm_exprs <- preprocessCore::normalize.quantiles(exprs_mx)

  colnames(norm_exprs) <- cnames
  rownames(norm_exprs) <- rnames

  norm_exprs
}


#' Normalize matrix
#'
#' @param exprs_dt data.table with raw, background-corrected expression. One
#' column per sample, one row per feature.
#' @param meta_data list of study-specific meta data
#'
#' @export
normalize_matrix <- function(exprs_dt, platform, verbose = FALSE) {

  badRows <- rowSums(is.na(exprs_dt)) > 0
  if ( sum(badRows) > 0 ) {
    warning("Removing ", sum(badRows), " rows with missing values")
    exprs_dt <- exprs_dt[complete.cases(exprs_dt)]
  }

  ## Prepare matrix

  feature_ids <- exprs_dt[ , feature_id ]
  exprs_dt[ , feature_id := NULL ]

  # Must ensure numeric values as conversion back to character can happen with
  # casting as matrix.
  exprs_mx <- as.matrix(apply(exprs_dt, 2, as.numeric))


  ## Call the correct method

  norm_exprs <- switch(
    platform,
    "Illumina" = normalize_microarray(exprs_mx,
                                      log2_transform = TRUE,
                                      verbose = verbose),
    "Affymetrix" = normalize_microarray(exprs_mx,
                                        log2_transform = FALSE,
                                        verbose = verbose),
    "NA" = normalize_rnaseq(exprs_mx,
                            verbose = verbose),
    "Stanford Functional Genomics Facility" = normalize_microarray(exprs_mx,
                                                                   log2_transform = FALSE,
                                                                   verbose = verbose),
    stop("Did not recognize platform: ", platform)
  )

  ## Convert back into data.table
  norm_exprs <- data.table(norm_exprs)
  norm_exprs[ , feature_id := feature_ids ]

  sample_id_columns <- grep("feature_id", names(norm_exprs), invert = TRUE)
  feature_id_column <- grep("feature_id", names(norm_exprs))
  setcolorder(norm_exprs, c(feature_id_column, sample_id_columns))

  norm_exprs
}

