suppressPackageStartupMessages({
  library(ImmuneSpaceR)
  library(HIPCMatrix)
  library(data.table)
  library(Rlabkey)
  library(testthat)
})


# Declare global test variables ------------------------------------------------

# Move netrc file out the way to test against login and password in .Renviron
if (!any(file.exists("~/.netrc", "~/_netrc"))) {
  assign("labkey.netrc.file", ImmuneSpaceR:::.get_env_netrc(), .GlobalEnv)
}

# Regardless of netrc presence, test against `ISR_machine` in .Renviron
assign("labkey.url.base", ImmuneSpaceR:::.get_env_url(), .GlobalEnv)

con_all <- CreateConnection("")
