
#' Process RNA-seq
#'
#' Somewhat redundant in case of GEO files as there is no actual processing
#' of raw counts files ...
#'
#' @param input_files input file names
#' @param study study accession eg \code{SDY269}
.processRnaSeq <- function(input_files, study){
  lf <- lapply(input_files, fread)
  exprs <- data.table(Reduce(f = function(x, y) {merge(x, y, all = TRUE)}, lf))
}
