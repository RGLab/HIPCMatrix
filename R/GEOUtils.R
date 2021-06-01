
#' Download Supplemental files from GEO
#'
#' Download, unzip, and return correct paths for GEO supplementary files
#'
#' @param accList accList
#' @param baseDir path to base directory (via \code{.getBaseDir()})
#' @param study study accession eg \code{SDY269}
#'
.dlSuppFls <- function(accList, baseDir, study){
  tmp <- sapply(accList, GEOquery::getGEOSuppFiles, makeDirectory = FALSE, baseDir = baseDir)
  fls <- list.files(baseDir)
  targetFlTerms <- "non-normalized|corrected|raw|cel|pbmc|count"
  rawFls <- fls[ grep(targetFlTerms, fls, ignore.case = TRUE) ]

  # Unzip any files if necessary - set `overwrite = TRUE` in case of processing fail
  flPaths <- rawFls[ grep("gz", rawFls) ]
  if (length(flPaths) > 0) {
    tmp <- sapply(file.path(baseDir, flPaths),
                  GEOquery::gunzip, overwrite = TRUE, remove = TRUE)
  }

  # find correct unzipped files with full paths
  # Do not include attempted intermediate file that may have been created if
  # run failed at a later point.
  fls <- file.path(baseDir, list.files(baseDir))
  fls <- fls[ grep(targetFlTerms, fls, ignore.case = TRUE)]
  inputFiles <- fls[ grep("gz|tar|RData", fls, invert = TRUE) ]
  inputFiles <- inputFiles[ grep(paste0(study, "_raw_expression"), inputFiles, invert = TRUE) ]
}

#' make id to GSM map
#'
#' Create map of GSM accessions to study-given IDs from GEO
#'
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
#' @param metaData list of study-specific meta data
#' @param study study accession eg \code{SDY269}
.makeIdToGsmMap <- function(gef, metaData, study){
  mp <- lapply(gef$geo_accession, function(x){
    tmp <- GEOquery::getGEO(x)
    nm <- tmp@header[[metaData$studyIdTerm]][[metaData$gsmMapIndex]]
    if (study %in% names(metaData$smplGsubTerms)) {
      nm <- gsub(metaData$smplGsubTerms[[study]]$old,
                 metaData$smplGsubTerms[[study]]$new,
                 nm)
    }
    return(c(x, nm))
  })
  mp <- data.frame(do.call(rbind, mp), stringsAsFactors = FALSE)
  colnames(mp) <- c("gsm","id")
  rownames(mp) <- NULL
  return(mp)
}

#' get GSM Supplementary Files
#'
#' Download Gsm supplementary files, unzip and return path
#'
#' @param gsm GEO GSM accession
#' @param baseDir path to base directory (via \code{.getBaseDir()})
.getGsmSuppFiles <- function(gsm, baseDir){
  info <- GEOquery::getGEOSuppFiles(gsm, makeDirectory = FALSE, baseDir = baseDir)
  GEOquery::gunzip(rownames(info), overwrite = TRUE, remove = TRUE)
  return( gsub("\\.gz", "", rownames(info)) )
}


#' Prep GEO files
#'
#' Generate flat files that are ready for processing from GEO "raw" data.
#' Warning! Raw data is highly variable for gse supplementary files
#'
#' @param study study accession eg \code{SDY269}
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
#' @param metaData list of study-specific meta data
#' @param inputFiles input file names
#'
#' @return path to raw, prepped input files
#'
.prepGeoFls <- function(study, gef, metaData, inputFiles){
  baseDir <- .getBaseDir(study, gef)

  # Case 1: raw data is in object returned by getGEO(gsm)
  # Only Illumina as of DR28
  if (metaData$dataInGsm == TRUE) {
    mxList <- lapply(gef$geo_accession, function(x){
      obj <- GEOquery::getGEO(x)
      tbl <- obj@dataTable@table
      tbl <- tbl[ , colnames(tbl) %in% c("ID_REF", metaData$gsmTblVarNm[[study]]) ]
      colnames(tbl)[[2]] <- x
      return(tbl)
    })
    inputFiles <- .mxListToFlatFile(mxList, baseDir, study)

  } else {
    if (metaData$useGsmSuppFls == TRUE) {
      # Case 2: raw data in gsm supp files - Illumina
      if (metaData$platform == "Illumina") {
        mxList <- lapply(gef$geo_accession, function(gsm){
          path <- .getGsmSuppFiles(gsm, baseDir)

          if (study %in% names(metaData$illuminaManifestFile)) {
            res <- limma::read.idat(idatfiles = path,
                                    bgxfile = metaData$illuminaManifestFile[[study]])
            em <- res$E
            pvals <- limma::detectionPValues(res)
            em <- data.table(gsm = em[,1], pvals = pvals[,1], ID_REF = res$genes$Probe_Id)
            em <- em[ !duplicated(em$ID_REF) ] # dups b/c single probe assigned to multiple array_ids
            setnames(em, "gsm", paste0(gsm, ".AVG_Signal"))
            setnames(em, "pvals", paste0(gsm, ".Detection Pval"))
          } else {
            em <- fread(path)
            em <- .subsetIlluminaEM(em)
            em <- .prepIlluminaHeaders(em)
            smplFormats <- "\\d{10}_[A-Z]"
            smplId <- regmatches(colnames(em)[[2]], regexpr(smplFormats, colnames(em)[[2]]))
            colnames(em) <- gsub(smplId, gsm, colnames(em))
          }
          return(em)
        })

        # Case 3: raw data in gsm supp files - Affymetrix
      } else if (metaData$platform == "Affymetrix") {
        inputFiles <- .dlSuppFls(accList = gef$geo_accession, baseDir, study)

        # Case 4: raw data in gsm supp files - Stanford custom HEEBO
      } else if (grepl("Stanford", metaData$platform)) {
        mxList <- lapply(gef$geo_accession, function(gsm){
          path <- .getGsmSuppFiles(gsm, baseDir)
          # Because of two colors, do background correction and processing here
          # to generate single expression value per probe
          em <- .processTwoColorArray(path)
          setnames(em, "gsm", gsm)
          return(em)
        })

        # Case 5: raw data in gsm supp files - RNAseq
      } else if (metaData$platform == "NA"){
        mxList <- lapply(gef$geo_accession, function(gsm){
          path <- .getGsmSuppFiles(gsm, baseDir)
          em <- fread(path)
          setnames(em, "V2", gsm) # ensemblId as 'V1'
          return(em)
        })
      }

      if (metaData$platform != "Affymetrix"){
        inputFiles <- .mxListToFlatFile(mxList, baseDir, study)
      }

    } else {
      # Cases 6 and 7: raw data in gse supp files - Illumina / RNAseq
      # temp handling for specialCase
      if(!metaData$useCustomRawFile){
        accList <- unique(unlist(lapply(gef$geo_accession, function(x){
          gsm <- getGEO(x)
          gse <- gsm@header$series_id
        })))
        inputFiles <- .dlSuppFls(accList, baseDir, study)
      }
      mxList <- lapply(inputFiles, fread)
      mxList <- .fixHeaders(mxList, study)

      # Case 6: Illumina raw data in gse supp files
      # Because multiple raw files need to be combined, must
      # address header issues prior to combination otherwise
      # untreated "Detection Pval" cols will cause dup error
      # during merge. Note: SDY400 handled in fixHeaders
      if (metaData$platform == "Illumina") {
        needMap <- study %in% names(metaData$smplGsubTerms) | metaData$gseNeedsMap
        if (needMap == TRUE) {
          mp <- .makeIdToGsmMap(gef, metaData, study)
        }

        mxList <- lapply(mxList, function(em){
          em <- .subsetIlluminaEM(em)
          em <- .prepIlluminaHeaders(em)
          if (needMap == TRUE) {
            # Fixed b/c some ids have escape char (e.g. ".")
            # Paste0 with "^" b/c some ids are numeric and confusable (e.g. "2.1" and "12.1")
            # lookahead to ensure full id before sep and not partial (e.g "PBMC_1" and "PBMC_12")
            # perl = TRUE for lookahead
            for (i in 1:nrow(mp)) {
              fixedId <- paste0("^", gsub(".", "\\.", mp$id[[i]], fixed = T), "(?=(\\.|_|$))")
              colnames(em) <- gsub(fixedId, mp$gsm[[i]], colnames(em), perl = TRUE)
            }
          }
          em <- em[ , grep("GSM|ID_REF", colnames(em)), with = FALSE]
        })
      }

      em <- Reduce(f = function(x, y) {merge(x, y)}, mxList)

      # Case 7: RNAseq in gse supp files
      # Header mapping assumes that names are in getGEO(gsm) object.
      # Need to check on a per study basis and tweak if need be.
      if (metaData$platform == "NA") {
        mp <- .makeIdToGsmMap(gef, metaData, study)
        em <- em[ , colnames(em) %in% c("GENES","V1", mp$id), with = FALSE ]
        nms <- grep("GENES|V1", colnames(em), invert = TRUE, value = TRUE)
        gsm <- mp$gsm[ match(nms, mp$id) ]
        setnames(em, nms, gsm)
      }

      inputFiles <- .writeSingleMx(em, baseDir, study)
    }
  }
  return(inputFiles)
}

