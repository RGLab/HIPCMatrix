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
#' \code{feature_id} column, with additional columns for each sample. Samples
#' should be labelled by biosample id.
#' @param feature_gene_map data.table mapping original gene id (eg probe id, ensembl id, gene alias)
#' to gene symbol. It should have two columns: \code{featureid} (original) and \code{genesymbol} (updated).
#' @param method Summarization method. \code{mean} or \code{max}.
#' \code{mean} will return the mean of all features mapping
#' to a gene alias for each sample. This is the method used in the ImmuneSpace
#' HIPCMatrix pipeline.
#' \code{max} will chose the probe which has the greatest geometric mean expression
#' across all samples. This is the method used in the ImmuneSignatures analysis.
#' @param verbose write verbose logging statements?
#'
#' @export
summarize_by_gene_symbol <- function(ge_dt,
                                     feature_gene_map,
                                     method = "mean",
                                     verbose = FALSE) {
  if (verbose) log_message("Summarizing matrix by gene symbol")
  em <- copy(ge_dt)
  # map feature to gene symbol
  em[, gene_symbol := feature_gene_map[match(em$feature_id, feature_gene_map$featureid), genesymbol]]
  # remove NA
  em <- em[!is.na(gene_symbol) & gene_symbol != "NA"]
  # Summarize
  if (method == "mean") {
    sum_exprs <- em[, lapply(.SD, mean),
      by = "gene_symbol",
      .SDcols = grep("^BS", colnames(em))
    ]
  } else if (method == "max") {
    em[, prb_avg := rowMeans(em[, grep("^BS", colnames(em)), with = FALSE])]
    maxPrb <- em[, prb_max := max(prb_avg), by = "gene_symbol"][prb_avg == prb_max]


    # Check and remove summary level duplicates w
    # here multiple probes have exact same average value AND same sample values
    maxPrb <- maxPrb[!(duplicated(maxPrb[, grep("feature_id", names(maxPrb), invert = TRUE), with = FALSE])), ]

    # handle cases where duplicated gs due to same average value,
    # but slightly different sample values. Keep last
    # e.g. SDY1289 - Lausanne Adult Cohort
    if (any(duplicated(maxPrb$gene_symbol))) {
      dup_gs <- maxPrb$gene_symbol[duplicated(maxPrb$gene_symbol)]
      for (dup in unique(dup_gs)) {
        rows <- grep(dup, maxPrb$gene_symbol)
        rm_rows <- rows[1:length(rows) - 1]
        maxPrb <- maxPrb[-rm_rows]
      }
    }

    sum_exprs <- maxPrb[, c("gene_symbol", grep("^BS", names(maxPrb), value = TRUE)), with = FALSE]
  }

  if (verbose) log_message(dim(ge_dt)[1], " features summarized to ", dim(sum_exprs)[1], " gene symbols")
  sum_exprs
}
