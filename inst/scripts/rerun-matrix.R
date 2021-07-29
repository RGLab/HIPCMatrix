# re-run create matrix workflow for an existing matrix
# already on the server.

# ----- Parse inputs -----
library(optparse)
option_list <- list(
  make_option(c("-m", "--matrix-name"),
              help = "Name of matrix stored in ImmuneSpace"),
  make_option(c("-l", "--labkey-url"),
              help = "labkey.url.base",
              default = "https://test.immunespace.org"),
  make_option(c("-r", "--reload"),
              help = "Force re-download of GEO files",
              action = "store_true",
              default = FALSE)
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
message(paste0(capture.output(opt), collaps = "\n"))

library(HIPCMatrix)
library(ImmuneSpaceR)
labkey.url.base <- opt$`labkey-url`
con <- CreateConnection("")
params <- HIPCMatrix:::get_input_params(con = con,
                                        matrix_name = opt$`matrix-name`)
params$verbose <- params$snapshot <- TRUE
params$reload <- opt$reload
do.call(runCreateMx,
        params)
