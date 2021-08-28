
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
  # Manually set curl options to ensure
  # that labkey.selectRows uses correct credential
  labkey.setCurlOptions(
    NETRC_FILE = ImmuneSpaceR:::.get_env_netrc()
  )
}

# Regardless of netrc presence, test against `ISR_machine` in .Renviron
assign("labkey.url.base", ImmuneSpaceR:::.get_env_url(), .GlobalEnv)

con_all <- HMX$new("")
SDY269 <- HMX$new("SDY269")
