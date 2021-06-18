#' Process RNA-seq
#'
#' Somewhat redundant in case of GEO files as there is no actual processing
#' of raw counts files ...
#'
#' @param input_files input file names
.process_rna_seq <- function(input_files){
  lf <- lapply(input_files, fread)
  exprs <- data.table(Reduce(f = function(x, y) {merge(x, y, all = TRUE)}, lf))
}
