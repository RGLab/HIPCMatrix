#' Function to migrate a copy of EMs from /Studies/ to a virtual study container
#'
#' @param ISserver immunespace server, either test or prod
#' @param virtualSdy subdirectory within HIPC container
#' @import Rlabkey
#' @export
#'

# Main Method
copyMatricesToVirtualSdy <- function(ISserver, virtualSdy){
  # Note: This should be run as the immunespace unix user on the RServe EC2 instance
  
  baseUrl <- ifelse( ISserver == "prod",
                     "https://www.immunespace.org",
                     "https://test.immunespace.org")

  runs <- labkey.selectRows(baseUrl = baseUrl,
                            folderPath = "/Studies/",
                            schemaName = "assay.ExpressionMatrix.matrix",
                            queryName = "Runs",
                            showHidden = TRUE)

  # subset to those in vSdy
  vSdyList <- list()
  vSdyList[["IS2"]] <- c("SDY56, SDY61, SDY63, SDY67, SDY80, SDY180, SDY212, SDY224, SDY269, SDY270, SDY400, SDY404, SDY520, SDY640, SDY984, SDY1119, SDY1260, SDY1264, SDY1276, SDY1289, SDY1291, SDY1293, SDY1294, SDY1325, SDY1328, SDY1364, SDY1368, SDY1370, SDY1373, SDY1529")
  vSdyList[["IS1"]]  <- c("SDY63, SDY67, SDY80, SDY212, SDY400, SDY404")
  studies <- vSdyList[[virtualSdy]]
  studies <- strsplit(studies, ", ")[[1]]
  sdy67 <- runs[ runs$Study == "SDY67", ]
  sdy67 <- sdy67[ grep("(B|b)atch", sdy67$Name, invert = virtualSdy != "IS1"), ]
  runs <- runs[ runs$Study %in% studies & runs$Study != "SDY67", ]
  runs <- rbind(sdy67, runs)

  # for each run cp over the normalized tsv only // 08.09.18 - will change and be unnecessary
  root <- "/share/files/Studies/"
  mid <- "/@files/analysis/exprs_matrices/"

  vPath <- paste0("/share/files/HIPC/", virtualSdy, "/@files/analysis/exprs_matrices/")

  for(i in 1:nrow(runs)){
    x <- runs[i, ]

    print(x$Name)

    mxPathOnSdy <- paste0(root, x$Study, mid, x$Name, ".tsv")
    runDir <- paste0(vPath, "Run", x$`Row Id`)
    mxPathOnVirt <- paste0(runDir, "/", x$Name, ".tsv")

    if (!dir.exists(runDir) ) {
      dir.create(runDir)
    }

    # the latest / normalized .tsv are not same?
    file.copy(from = mxPathOnSdy,
              to = mxPathOnVirt,
              overwrite = TRUE)

    # if IS1
    if(virtualSdy == "IS1"){
      file.copy(from = paste0(mxPathOnSdy, ".immsig"),
                to = paste0(mxPathOnVirt, ".immsig"),
                overwrite = FALSE)
    }

    file.copy(from = paste0(mxPathOnSdy, ".raw"),
              to = paste0(mxPathOnVirt, ".raw"),
              overwrite = TRUE)

    file.copy(from = paste0(mxPathOnSdy, ".summary"),
              to = paste0(mxPathOnVirt, ".summary"),
              overwrite = TRUE)

    file.copy(from = paste0(mxPathOnSdy, ".summary.orig"),
              to = paste0(mxPathOnVirt, ".summary.orig"),
              overwrite = TRUE)

    # # RM any files without .tsv (due to deprecated create-matrix.xml)
    # fls <- list.files(runDir)
    # rmFls <- fls[ grep("tsv|xml", fls, invert = T)]
    # rmFls <- file.path(runDir, rmFls)
    # done <- sapply(rmFls, file.remove)
  }
}


