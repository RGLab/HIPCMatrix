#' Function to migrate a copy of FA / FAS from /Studies/ to a virtual study container
#'
#' @param ISserver immunespace server, either test or prod
#' @param virtualSdy subdirectory within HIPC container
#' @param fasGrep grep statement to use in subsetting FAS to add just a few new sets
#' @param verbose prints name of FAS being copied
#' @import Rlabkey
#' @export
#'

# Main Method
addAnnoToVirtualSdy <- function(ISserver, virtualSdy, fasGrep = NULL, verbose = FALSE){

    # Convert server to baseUrl
    baseUrl <- ifelse( ISserver == "prod",
                     "https://www.immunespace.org",
                     "https://test.immunespace.org")

    # Path assumed to be from HIPC directory
    vSdyPath <- paste0("/HIPC/", virtualSdy)

    # retrieve FAS
    fas <- labkey.selectRows(baseUrl = baseUrl,
                             folderPath = "/Studies/",
                             schemaName = "microarray",
                             queryName = "FeatureAnnotationSet",
                             colNameOpt = "fieldname",
                             showHidden = TRUE)
    if (!is.null(fasGrep)) {
      fas <- fas[ grep(fasGrep, fas$Name), ]
    }
    fas <- fas[ order(fas$RowId), ]
    toImport <- data.frame(fas, stringsAsFactors = FALSE)
    toImport[is.na(toImport)] <- ""
    toImport <- toImport[ grep("RowId", colnames(toImport), invert = TRUE)]

    # Get container ID
    container <- labkey.selectRows(baseUrl = baseUrl,
                                    folderPath = vSdyPath,
                                    schemaName = "core",
                                    queryName = "containers",
                                    colNameOpt = "fieldname",
                                    colSelect = "EntityId")
    toImport$Container <- container[1,1]

    # Import FAS to vSdy
    impFas <- labkey.importRows(baseUrl = baseUrl,
                                folderPath = vSdyPath,
                                schemaName = "microarray",
                                queryName = "FeatureAnnotationSet",
                                toImport = toImport)

    # Get new rowId mappings
    newFas <- labkey.selectRows(baseUrl = baseUrl,
                                folderPath = vSdyPath,
                                schemaName = "microarray",
                                queryName = "FeatureAnnotationSet",
                                colNameOpt = "fieldname",
                                showHidden = TRUE)
    if (!is.null(fasGrep)) {
      newFas <- newFas[ grep(fasGrep, newFas$Name), ]
    }
    newFas <- newFas[ order(match(newFas$Name, fas$Name)), ]

    if (!all(fas$Name == newFas$Name)) {
      stop("Imported FAS is not the same as /Studies/. Please fix.")
    }

    rowMap <- data.frame(old = fas$RowId,
                         new = newFas$RowId,
                         stringsAsFactors = FALSE)

    # Update the featureAnnotation and Import
    # NOTE: colSelect "all" option aka "*" used in order to get FASid
    fasIdFilt <- makeFilter(c('FeatureAnnotationSetId', 'IN', paste(fas$RowId, collapse = ';')))
    fa <- labkey.selectRows(baseUrl = baseUrl,
                            folderPath = "/Studies/",
                            schemaName = "microarray",
                            queryName = "FeatureAnnotation",
                            colNameOpt = "fieldname",
                            colSelect = "*",
                            colFilter = fasIdFilt,
                            showHidden = TRUE)
    newFa <- fa[ , colnames(fa) %in% c("Container","FeatureAnnotationSetId", "FeatureId", "GeneSymbol")]
    newFa$Container <- unique(toImport$Container)
    newFa$FeatureAnnotationSetId <- rowMap$new[ match(newFa$FeatureAnnotationSetId, rowMap$old)]
    newFa[is.na(newFa)] <- ""
    
    # Importing > 100k rows bogs down, so do subsets
    numCuts <- ceiling(nrow(newFa)/100000)
    if(verbose){
      print(paste("Total Rows to Import:", nrow(newFa)))
    }
    
    for(i in 1:numCuts){
      low <- 100000 * (i-1) + 1
      high <- 100000 * i
      if(high > nrow(newFa)){
        high <- nrow(newFa)
      }
      
      if(verbose){
        print(paste0("low: ", low, " high: ", high))
      }
      
      tmp <- newFa[low:high,]
      doneFa <- labkey.importRows(baseUrl = baseUrl,
                                  folderPath = vSdyPath,
                                  schemaName = "microarray",
                                  queryName = "FeatureAnnotation",
                                  toImport = tmp)
    }
}
