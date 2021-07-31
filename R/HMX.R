#' HIPCMatrix Connection
#'
#' @description ImmuneSpace connection object with additional methods
#' for gene expression analysis
#'
#' @importFrom R6 R6Class
#' @export
HMX <- R6::R6Class(
  classname = "HIPCMatrixConn",
  inherit = ImmuneSpaceR:::ISCon
)
