getMetaData_old <- function(study, gef, fasId = NULL, con = NULL) {
  if (is.null(fasId)) {
    fasId <- unique(
      Rlabkey::labkey.selectRows(baseUrl = con$config$labkey.url.base,
                                      folderPath = con$config$labkey.url.path,
                                      schemaName = "assay.ExpressionMatrix.matrix",
                                      queryName = "Runs",
                                      colNameOpt = "fieldname")$featureSet
    )
  }
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
  metaData$specialCase <- study == "SDY1529" & (0 %in% unique(gef$study_time_collected))
  metaData$useCustomRawFile <- study %in% c("SDY224", "SDY1324") | metaData$specialCase

  metaData
}



# ----- spec for return from getMetaData_new: -----
# list(
#   raw_file_location = "gsm|gse|gsm_supp_files|immport|custom",
#   id_to_gse_mapping_info = NULL | list(
#     study_id_term = "description|title",
#     gsm_map_index = 1|2,
#     id_regex_map = NULL | list(
#       old = "regex",
#       new = "regex"
#     )
#   ),
#   platform = "affymetrix|illumina|RNAseq",
#   illumna_manifest_file = NULL | "filepath",
#   gsm_table_var_name = NULL | "Value"
# )


#' Get Meta Data
#'
#' @param study study accession
#' @param gef gene expression files
#' @param fasId Feature Annotation Set ID (from FeatureAnnotationSet table)
#'
#' @export
getMetaData <- function(study, gef, fasId) {

  # Locations for files for all studies.
  # There is one study: SDY1529, which has a custom file for one matrix.
  # All other studies use the same location for all matrices.

  # NOTES: Preference is to use GEO when selectedBiosamples have both gsm and
  # flat files because there is often info in GEO about whether files are raw
  # or normalized that is not available in ImmPort. However, there are exceptions
  # and those are handled here - e.g. SDY212 where we know from ImmSig1 that the
  # available flat file is truly non-normalized data.
  file_locations <- list(

    # raw data is kept in the eset provided by getGEO()
    gsm_soft = "SDY1289",

    # gsm_supp_files: Usually refers to Affymetrix studies that have the CEL.gz files
    # loaded into GEO as a supplementary file to the single GSM accession as opposed to
    # the Illumina that often have a single tsv.gz file in the series accession (GSE)
    gsm_supp_files = c("SDY80", "SDY113", "SDY180", "SDY269",
                       "SDY406", "SDY984", "SDY1260", "SDY1264",
                       "SDY1293", "SDY270", "SDY1291", "SDY212",
                       "SDY315", "SDY305", "SDY1328", "SDY1368",
                       "SDY1370", "SDY1119", "SDY1294", "SDY1256",
                       "SDY1412", "SDY1267", "SDY1086"),

    # gse_supp_files: Mostly illumina and RNA-seq. All data is included in
    # GSE supp files.
    gse_supp_files = c("SDY61", "SDY144", "SDY63", "SDY400", "SDY404",
                       "SDY520", "SDY640", "SDY56", "SDY789", "SDY820",
                       "SDY522", "SDY888", "SDY1276",
                       "SDY1325", "SDY1364", "SDY1373", "SDY1092",
                       "SDY903", "SDY787", "SDY1529"),

    # Use raw files from ImmPort. Be careful because these may already be
    # normalized.
    immport = c("SDY1630", "SDY645", "SDY376", "SDY299",
                "SDY89", "SDY67", "SDY28", "SDY34",
                "SDY112", "SDY690", "SDY597", "SDY387",
                "SDY372", "SDY368", "SDY364", "SDY312", "SDY301",
                "SDY296", "SDY667", "SDY300", "SDY162"),

    # custom: Use some other custom file. Will have to be pulled from
    # ImmuneSpace server. This is generally a custom file provided by
    # immport as they wait for files to be updated in GEO or immport
    custom = c("SDY224", "SDY1324", "SDY1529"),

    # No raw data in any location
    not_available = c()
  )


  # Info for mapping ids found in GSE file to individual biosample ids (GSM)
  id_to_gsm_mapping_info <- list(
    gsm_needs_map = c("SDY404", "SDY522", "SDY1325", "SDY1364", "SDY144",
                       "SDY400", "SDY640", "SDY520", "SDY1529", "SDY1276",
                       "SDY224", "SDY63", "SDY888", "SDY1373", "SDY180",
                       "SDY787"),

    # Term from getGEO(gsm) header to use for extracting sample id
    use_gsm_description = c("SDY144", "SDY180", "SDY522", "SDY1373", "SDY1364",
                            "SDY1325", "SDY640", "SDY520"),

    # Index of samplename in vector from getGEO()
    use_gsm_index_2 = c("SDY180", "SDY640", "SDY520"),

    # Gsub terms for mapping id to gse
    id_regex_map_list = list(
      SDY1276 = list(old = "WholeBloodRNA_", new = ""),
      SDY224  = list(old = " \\[PBMC\\]", new = ""),
      SDY63   = list(old = "^101", new = "10"),
      SDY888  = list(old = "( |)_((N|n)egative|(S|s)econdary)", new = "_RNASeq"),
      SDY1373 = list(old = "Sample name: ", new = ""),
      SDY180  = list(old = "([0-9])([A-Z])", new = "\\1_\\2"),
      SDY787  = list(old = "\\D", new = "") # Replace all non-digits
    )
  )

  #**illuminaManifestFile**: for studies with Illumina idat files that need bgx
  # manifest files.  These are found through the Illumina website and stored in
  # the UpdateAnno package. Below creates a temp file to store this data.
  # TODO: assign file based on fasId
  illumina_manifest_files <- list(
    SDY1368 = "HumanHT-12_V4_0_R2_15002873_B.bgx"
  )


  # Custom list of raw values column name for gsm-based data
  gsm_table_var_name <- list(
    SDY1289 = "AVERAGE_SIGNAL",
    SDY1293 = "VALUE"
  )

  # ----- Create custom metaData list for study -----

  metaData <- list()


  # Set File location
  if ( study == "SDY1529" ) {
    metaData$file_location <- ifelse( (0 %in% unique(gef$study_time_collected)),
                                      "custom",
                                      "gse_supp_files")
  } else {
    file_location <- names(file_locations)[
      vapply(file_locations, function(location) study %in% location,
             FUN.VALUE = TRUE)
    ]
    if (length(file_location) != 1) {
      stop("Could not determine location of input files.")
    }
    metaData$file_location <- file_location
  }

  # Set mapping info (Will be NULL if no mapping needed)
  if ( study %in% id_to_gsm_mapping_info$gsm_needs_map ) {
    metaData$id_to_gse_mapping_info <- list(
      study_id_term = ifelse(study %in% id_to_gsm_mapping_info$use_gsm_description, "description", "title"),
      gsm_map_index = ifelse(study %in% id_to_gsm_mapping_info$sue_gsm_index_2, 2, 1),
      id_regex_map  = id_to_gsm_mapping_info$id_regex_map_list[[study]]
    )
  }

  # Set platform (derive from fasId)

  # **platform**: sequencing platform (Affymetrix, Illumina, or 'NA', aka RNAseq)
  fas <- data.table(labkey.selectRows(baseUrl = labkey.url.base,
                                      folderPath = "/Studies/",
                                      schemaName = "Microarray",
                                      queryName = "FeatureAnnotationSet",
                                      colNameOpt = "rname",
                                      showHidden = T))
  metaData$platform <- fas$vendor[ fas$rowid == fasId ]

  metaData$illumina_manifest_file <- illumina_manifest_files[[study]]

  metaData$gsm_table_var_name <- gsm_table_var_name[[study]]

  metaData
}

