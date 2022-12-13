# re-run create matrix workflow for all existing matrices
# on the server

# ----- Parse inputs -----
library(optparse)
option_list <- list(
  make_option(c("-l", "--labkey-url"),
              help = "labkey.url.base",
              default = "https://datatools-dev.immunespace.org")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
message(paste0(capture.output(opt), collaps = "\n"))

library(HIPCMatrix)
library(ImmuneSpaceR)
labkey.url.base <- opt$`labkey-url`
con <- CreateConnection("")
matrices <- con$listGEMatrices()$name
res <- lapply(matrices,

              function(matrix_name) {
                HIPCMatrix:::log_message("----- ", matrix_name, " -----")
                result <- try({
                  params <- HIPCMatrix:::get_input_params(con = con,
                                                          matrix_name = opt$`matrix-name`)
                  params$verbose <- params$snapshot <- TRUE
                  do.call(runCreateMx,
                          params)
                }, silent = TRUE)
                if ("try-error" %in% class(result)) {
                  HIPCMatrix:::log_message(matrix_name, " failed: ", result)
                  return(res)
                } else {
                  HIPCMatrix:::log_message("Success!")
                  return(TRUE)
                }
              }
)

failures <- matrices[!isTRUE(unlist(res))]
if (length(failures) > 0) {
  result_log <- paste0("/share/files/Studies/@files/rerun_all_matrices_",
                       strftime(Sys.Date(), "%Y%m%d"),
                       ".rds")
  writeRDS(res, result_log)
  HIPCMatrix:::log_message(length(failures), " matrices failed:")
  lapply(failures, HIPCMatrix:::log_message)
  HIPCMatrix:::log_message("See ", result_log, " for details.")
} else {
  HIPCMatrix:::log_message("All matrices successfully processed! Horray!")
}
