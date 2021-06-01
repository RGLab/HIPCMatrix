
#' Process RNA-seq
#'
#' Somewhat redundant in case of GEO files as there is no actual processing
#' of raw counts files ...
#'
#' @param inputFiles input file names
#' @param study study accession eg \code{SDY269}
.processRnaSeq <- function(inputFiles, study){
  lf <- lapply(inputFiles, fread)
  exprs <- data.table(Reduce(f = function(x, y) {merge(x, y, all = TRUE)}, lf))
}
