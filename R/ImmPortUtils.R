#' Prep ImmPort Files
#'
#' @param study study accession eg \code{SDY269}
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
#' @param meta_data list of study-specific meta data
#' @param input_files input file names
#'
#' @return path to raw, prepped input files
#'
.prep_immport_files <- function(study,
                                gef,
                                platform,
                                input_files,
                                analysis_dir){

  supp_files_dir <- .get_supp_files_dir(analysis_dir,
                                        gef)

  if( platform == "Illumina" ) {
    ge_list <- lapply(input_files, function(path){
      raw_illumina_dt <- fread(path)
      raw_illumina_dt <- .subset_raw_illumina_dt(raw_illumina_dt)
      raw_illumina_dt <- .prep_illumina_headers(raw_illumina_dt)
    })
  } else
    if ( platform == "NA" ) {
      ge_list <- lapply(input_files, fread)
      ge_list <- .fix_headers(ge_list, study)
    }

  input_files <- .ge_list_to_flat_file(ge_list, supp_files_dir, study)
  return(input_files)
}



#' Map experiment-sample or geo accessions to biosample accessions
#'
#' @param exprs_dt data.table of gene expression with one column per sample,
#' one row per feature.
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
#'
.map_sampleid_to_biosample_accession <- function(exprs_dt, gef){
  if ( any(grepl("^(ES|GSM)\\d{6,7}$", colnames(exprs_dt))) ) {
    col_to_use <- ifelse( any(grepl("ES", colnames(exprs_dt))),
                        "expsample_accession",
                        "geo_accession")
  } else {
    col_to_use <- "file_info_name"
  }
  sampleids <- grep("feature_id", colnames(exprs_dt), value = TRUE, invert = TRUE)
  biosample_accessions <- gef$biosample_accession[ match(sampleids, gef[[col_to_use]]) ]

  # remove samples without matching biosample accession
  sampleids <- sampleids[!is.na(biosample_accessions)]
  biosample_accessions <- biosample_accessions[!is.na(biosample_accessions)]

  setnames(exprs_dt, sampleids, biosample_accessions)
  return(exprs_dt)
}
