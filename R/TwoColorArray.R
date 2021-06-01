
#' Process Two Color Array
#'
#' Two color array processing using limma and assuming genepix files.
#' bc.method is the background correction method and normexp is used to match
#' work with Illumina and Affymetrix.
.processTwoColorArray <- function(path){
  RG <- limma::read.maimages(files = path, source = "genepix")
  MA <- limma::normalizeWithinArrays(RG, bc.method = "normexp", method = "none")
  em <- data.table(feature_id = MA$genes$ID, gsm = MA$A[,1])

  # RM dup probes due to multiple probes per spot
  # RM NA vals possibly from background correction issues
  em <- em[ !duplicated(feature_id) & !is.na(gsm) ]

  # RM unmappable feature_ids
  em <- em[ grep("EMPTY|bsid", feature_id, invert = TRUE) ]
}
