#' Get Meta Data
#'
#' @param study study accession
#' @param gef gene expression files
#' @param fasId Feature Annotation Set ID (from FeatureAnnotationSet table)
getMetaData <- function(study, gef, fasId) {
  # Manually curate list object to avoid hardcoding elsewhere
  metaData <- list()

  # **isGeo**: Preference is to use GEO when selectedBiosamples have both gsm and
  # flat files because there is often info in GEO about whether files are raw
  # or normalized that is not available in ImmPort. However, there are exceptions
  # and those are handled here - e.g. SDY212 where we know from ImmSig1 that the
  # available flat file is truly non-normalized data.
  dontUseGeo <- study %in% c("SDY224")
  metaData$isGeo <- any(!is.na(gef$geo_accession)) & !dontUseGeo

  # **dataInGsm**: means that the raw data is kept in the eset provided by getGEO()
  metaData$dataInGsm <- study %in% c("SDY1289")

  # **useGsmSuppFls**: Usually refers to Affymetrix studies that have the CEL.gz files
  # loaded into GEO as a supplementary file to the single GSM accession as opposed to
  # the Illumina that often have a single tsv.gz file in the series accession (GSE)
  metaData$useGsmSuppFls <- study %in% c("SDY80", "SDY113", "SDY180", "SDY269",
                                         "SDY406", "SDY984", "SDY1260", "SDY1264",
                                         "SDY1293", "SDY270", "SDY1291", "SDY212",
                                         "SDY315", "SDY305", "SDY1328", "SDY1368",
                                         "SDY1370", "SDY1119", "SDY1294", "SDY1256",
                                         "SDY1412", "SDY1267", "SDY1086")

  #**illuminaManifestFile**: for studies with Illumina idat files that need bgx
  # manifest files.  These are found through the Illumina website and stored in
  # the UpdateAnno package. Below creates a temp file to store this data.
  # TODO: assign file based on fasId
  metaData$illuminaManifestFile <- list(
    SDY1368 = "HumanHT-12_V4_0_R2_15002873_B.bgx"
  )

  if (study %in% names(metaData$illuminaManifestFile)){
    manifestUrl <- paste0("https://github.com/RGLab/UpdateAnno/raw/main/CreateMatrixAssets/IlluminaManifests/", metaData$illuminaManifestFile[[study]])
    tmpFl <- tempfile()
    download.file(url = manifestUrl, destfile = tmpFl, quiet = TRUE)
    metaData$illuminaManifestFile[[study]] <- tmpFl
  }

  # **studyIdTerm**: For extracting sample id from getGEO(gsm) object
  useDescription <- study %in% c("SDY144", "SDY180", "SDY522", "SDY1373", "SDY1364",
                                 "SDY1325", "SDY640", "SDY520")
  metaData$studyIdTerm <- ifelse(useDescription, "description", "title")

  # **smplGsubTerms**: Custom gsub terms for allowing the mapping of study-given ids
  # from the supplementary files to the ids found in GEO in the header object.
  metaData$smplGsubTerms <- list(
    SDY1276 = list(old = "WholeBloodRNA_", new = ""),
    SDY224  = list(old = " \\[PBMC\\]", new = ""),
    SDY63   = list(old = "^101", new = "10"),
    SDY888  = list(old = "( |)_((N|n)egative|(S|s)econdary)", new = "_RNASeq"),
    SDY1373 = list(old = "Sample name: ", new = ""),
    SDY180  = list(old = "([0-9])([A-Z])", new = "\\1_\\2"),
    SDY787  = list(old = "\\D", new = "") # Replace all non-digits
  )

  # **gseNeedsMap**: Studies that need id-to-gsm mapping from gse supp files
  # without special gsub terms
  metaData$gseNeedsMap <- study %in% c("SDY404", "SDY522", "SDY1325", "SDY1364", "SDY144",
                                       "SDY400", "SDY640", "SDY520", "SDY1529")

  # **gsmMapIndex**: Index of samplename in vector from getGEO()
  useSecond <- c("SDY180", "SDY640", "SDY520")
  metaData$gsmMapIndex <- ifelse( study %in% useSecond, 2, 1)

  # **gsmTblVarNm**: Custom list of raw values column name for gsm-based data
  metaData$gsmTblVarNm <- list(
    SDY1289 = "AVERAGE_SIGNAL",
    SDY1293 = "VALUE"
  )

  # **noRaw**: No raw data available in GEO, ImmPort, or customRawFile
  metaData$noRaw <- study %in% c()

  # **platform**: sequencing platform (Affymetrix, Illumina, or 'NA', aka RNAseq)
  fas <- data.table(labkey.selectRows(baseUrl = labkey.url.base,
                                      folderPath = "/Studies/",
                                      schemaName = "Microarray",
                                      queryName = "FeatureAnnotationSet",
                                      colNameOpt = "rname",
                                      showHidden = T))
  metaData$platform <- fas$vendor[ fas$rowid == fasId ]

  # **useCustomRawFile**: For some studies a custom file has been provided by ImmPort
  # temporarily while they update a study. These should be checked periodically.
  specialCase <- study == "SDY1529" & (0 %in% unique(gef$study_time_collected))
  metaData$useCustomRawFile <- study %in% c("SDY224", "SDY1324") | specialCase

  metaData
}
