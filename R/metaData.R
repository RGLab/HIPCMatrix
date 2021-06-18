#' List of study-specific metadata
#'
#'
#' @format A list with 4 entries:
#' \describe{
#'   \item{file_locations}{list of file locations, with a vector of studies
#'   for each of \code{gsm_soft}, \code{gsm_supp_files}, \code{gse_supp_files},
#'   \code{immport}, \code{custom}, \code{not_available}}
#'   \item{id_to_gsm_mapping_info}{list of mapping info for studies which have
#'   data in the GSE and need additional info to map back to the GSM per sample.
#'   Includes: \code{gsm_needs_map}, \code{use_gsm_description},
#'   \code{use_gsm_index_2}, \code{id_regex_map_list}}
#'   \item{illumina_manifest_files}{List with one entry per study which
#'   has an illumina manifest file}
#'   \item{gsm_table_var_name}{list of raw values column name for gsm-based data}
#' }
"meta_data_list"


# ----- format for return from get_meta_data: -----
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
#' @param fas_id Feature Annotation Set ID (from FeatureAnnotationSet table)
#'
#' @export
get_meta_data <- function(study, gef, fas_id) {

  # ----- Create custom meta_data list for study -----

  meta_data <- list()

  # Set File location
  if ( study == "SDY1529" ) {
    meta_data$file_location <- ifelse( (0 %in% unique(gef$study_time_collected)),
                                      "custom",
                                      "gse_supp_files")
  } else {
    file_location <- names(meta_data_list$file_locations)[
      vapply(meta_data_list$file_locations, function(location) study %in% location,
             FUN.VALUE = TRUE)
    ]
    if (length(file_location) != 1) {
      stop("Could not determine location of input files.")
    }
    meta_data$file_location <- file_location
  }

  # Set mapping info (Will be NULL if no mapping needed)
  if ( study %in% meta_data_list$id_to_gsm_mapping_info$gsm_needs_map ) {
    meta_data$id_to_gse_mapping_info <- list(
      study_id_term = ifelse(study %in% meta_data_list$id_to_gsm_mapping_info$use_gsm_description, "description", "title"),
      gsm_map_index = ifelse(study %in% meta_data_list$id_to_gsm_mapping_info$use_gsm_index_2, 2, 1),
      id_regex_map  = meta_data_list$id_to_gsm_mapping_info$id_regex_map_list[[study]]
    )
  }

  # Set platform (derive from fas_id)

  # **platform**: sequencing platform (Affymetrix, Illumina, or 'NA', aka RNAseq)
  fas <- data.table(labkey.selectRows(baseUrl = labkey.url.base,
                                      folderPath = "/Studies/",
                                      schemaName = "Microarray",
                                      queryName = "FeatureAnnotationSet",
                                      colNameOpt = "rname",
                                      showHidden = T))
  meta_data$platform <- fas$vendor[ fas$rowid == fas_id ]

  meta_data$illumina_manifest_file <- meta_data_list$illumina_manifest_files[[study]]

  meta_data$gsm_table_var_name <- meta_data_list$gsm_table_var_name[[study]]

  meta_data
}

