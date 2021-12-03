meta_data_list <- list(

  # Locations for files for all studies.
  # There is one study: SDY1529, which has a custom file for one matrix.
  # All other studies use the same location for all matrices.

  # NOTES: Preference is to use GEO when selected_biosamples have both gsm and
  # flat files because there is often info in GEO about whether files are raw
  # or normalized that is not available in ImmPort. However, there are exceptions
  # and those are handled here - e.g. SDY212 where we know from ImmSig1 that the
  # available flat file is truly non-normalized data.
  file_locations = list(

    # raw data is kept in the eset provided by getGEO()
    gsm_soft = "SDY1289",

    # gsm_supp_files: Usually refers to Affymetrix studies that have the CEL.gz files
    # loaded into GEO as a supplementary file to the single GSM accession as opposed to
    # the Illumina that often have a single tsv.gz file in the series accession (GSE)
    gsm_supp_files = c(
      "SDY56", "SDY61", "SDY80", "SDY113", "SDY180",
      "SDY269", "SDY1086", "SDY1412", "SDY1267",
      "SDY406", "SDY984", "SDY1260", "SDY1264",
      "SDY1293", "SDY270", "SDY1291", "SDY212",
      "SDY315", "SDY305", "SDY1328", "SDY1368",
      "SDY1370", "SDY1119", "SDY1294", "SDY1256"
    ),

    # gse_supp_files: Mostly illumina and RNA-seq. All data is included in
    # GSE supp files.
    gse_supp_files = c(
      "SDY144", "SDY63", "SDY400", "SDY404",
      "SDY520", "SDY640", "SDY789", "SDY820",
      "SDY522", "SDY888", "SDY1276",
      "SDY1325", "SDY1364", "SDY1373", "SDY1092",
      "SDY903", "SDY787", "SDY1529"
    ),

    # Use raw files from ImmPort. Be careful because these may already be
    # normalized.
    immport = c(
      "SDY1630", "SDY645", "SDY376", "SDY299",
      "SDY89", "SDY67", "SDY28", "SDY34",
      "SDY112", "SDY690", "SDY597", "SDY387",
      "SDY372", "SDY368", "SDY364", "SDY312", "SDY301",
      "SDY296", "SDY667", "SDY300", "SDY162"
    ),

    # custom: Use some other custom file. Will have to be pulled from
    # ImmuneSpace server. This is generally a custom file provided by
    # immport as they wait for files to be updated in GEO or immport.
    # SDY1529 has a custom file for baseline samples only.
    custom = c("SDY224", "SDY1324", "SDY1529"),

    # No raw data in any location
    not_available = c()
  ),
  custom_file_info = list(
    SDY224 = list(
      directory = "",
      file_identifier_regex = "reads"
    ),
    SDY1324 = list(
      directory = "raw_counts",
      file_identifier_regex = "RawCounts"
    ),
    SDY1529 = list(
      directory = "author_data",
      file_identifier_regex = "GA"
    )
  ),


  # Info for mapping ids found in GSE file to individual biosample ids (GSM)
  id_to_gsm_mapping_info = list(
    gsm_needs_map = c(
      "SDY404", "SDY522", "SDY1325", "SDY1364", "SDY144",
      "SDY400", "SDY640", "SDY520", "SDY1529", "SDY1276",
      "SDY224", "SDY63", "SDY888", "SDY1373", "SDY180",
      "SDY787"
    ),

    # Term from getGEO(gsm) header to use for extracting sample id
    use_gsm_description = c(
      "SDY144", "SDY180", "SDY522", "SDY1373", "SDY1364",
      "SDY1325", "SDY640", "SDY520"
    ),

    # Index of samplename in vector from getGEO()
    use_gsm_index_2 = c("SDY180", "SDY640", "SDY520"),

    # Gsub terms for mapping id to gse
    id_regex_map_list = list(
      SDY1276 = list(old = "WholeBloodRNA_", new = ""),
      SDY224 = list(old = " \\[PBMC\\]", new = ""),
      SDY63 = list(old = "^101", new = "10"),
      SDY888 = list(old = "( |)_((N|n)egative|(S|s)econdary)", new = "_RNASeq"),
      SDY1373 = list(old = "Sample name: ", new = ""),
      SDY180 = list(old = "([0-9])([A-Z])", new = "\\1_\\2"),
      SDY787 = list(old = "\\D", new = "") # Replace all non-digits
    )
  ),

  #** illuminaManifestFile**: for studies with Illumina idat files that need bgx
  # manifest files.  These are found through the Illumina website and stored in
  # the UpdateAnno package. Below creates a temp file to store this data.
  # TODO: assign file based on fas_id
  illumina_manifest_files = list(
    SDY1368 = "HumanHT-12_V4_0_R2_15002873_B.bgx"
  ),


  # Custom list of raw values column name for gsm-based data
  gsm_table_var_name = list(
    SDY1289 = "AVERAGE_SIGNAL",
    SDY1293 = "VALUE"
  )
)
usethis::use_data(meta_data_list, overwrite = TRUE)
