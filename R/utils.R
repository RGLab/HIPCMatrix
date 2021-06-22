#' Get supplemental files directory
#'
#' Create run-specific base directory for supplementary files.
#' Important that run specific directory is created to avoid overlap with previous runs.
#' This done with unique cell_type * arm_accession string (aka cohort_type)
#'
#' @param analysis_dir analysis directory
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
.get_supp_files_dir <- function(analysis_dir,
                                gef){
  supp_files_dir <- file.path(analysis_dir,
                              "supp_files",
                              paste0(unique(gef$type),
                                     "_",
                                     unique(gef$arm_accession)))
  dir.create(supp_files_dir, recursive = TRUE)
  return(supp_files_dir)
}

#' Write raw expression to a file
#'
#' Writes \code{ge_tbl} to a (tsv-formatted) \code{.txt} file, using controlled
#' format and file name.
#'
#' @return path to raw expression file
#'
#' @param ge_tbl object containing gene expression values (matrix, data.frame, etc)
#' @param supp_files_dir path to base directory (via \code{.get_supp_files_dir()})
#' @param study study accession eg \code{SDY269}
.write_raw_expression <- function(ge_tbl, supp_files_dir, study) {
  input_files <- file.path(supp_files_dir, paste0(study, "_raw_expression.txt"))
  fwrite(ge_tbl, file = input_files, sep = "\t", quote = FALSE, row.names = FALSE)
  return(input_files)
}

#' fix headers
#'
#' Custom fixes for individual studies
#'
#' @param ge_list list of matrices for one study
#' @param study study accession eg \code{SDY269}
.fix_headers <- function(ge_list, study){
  if (study == "SDY224") {
    ge_list <- lapply(ge_list, function(x){
      setnames(x, as.character(x[1,]))
      x <- x[-(1:2),]
      colnames(x)[[1]] <- "ID_REF"
      return(x)
    })
  } else if (study == "SDY400") {
    # Using mapping file provided by Hailong Meng at Yale, Dec 2018
    # since note in header of file is misleading due to gsm swaps made
    # later based on knowledge of switched samples.
    ge_list <- lapply(ge_list, function(x){
      mp <- fread("/share/files/Studies/SDY400/@files/rawdata/gene_expression/SDY400_HeaderMapping.csv")
      setnames(x, colnames(x), as.character(x[2,]))
      x <- x[-(1:2),]
      smpls <- grep("SAMPLE", colnames(x), value = T)
      titles <- mp$Title[ match(smpls, mp$Sample) ]
      setnames(x, smpls, titles)
      return(x)
    })
  } else if (study == "SDY1325") {
    ge_list <- lapply(ge_list, function(x){
      setnames(x, colnames(x), as.character(x[5,]))
      x <- x[6:nrow(x),]
      return(x)
    })
  } else if (study == "SDY1324") {
    # Custom header mapping provided by authors via P.Dunn Dec 2018.
    ge_list <- lapply(ge_list, function(x){
      mp <- fread("/share/files/Studies/SDY1324/@files/rawdata/gene_expression/raw_counts/SDY1324_Header_Mapping.csv")
      accs <- grep("V1", colnames(x), invert = TRUE, value = TRUE)
      esNms <- mp$experimentAccession[ match(accs, mp$AuthorGivenId) ]
      setnames(x, accs, esNms)
      return(x)
    })
  } else if (study == "SDY787") {
    # Fist number is unique id
    ge_list <- lapply(ge_list, function(x) {
      # Remove first "_" and everything following
      setnames(x, colnames(x), gsub("_.*$", "", colnames(x)))
    })
  }
  return(ge_list)
}

#' matrix list to flat file
#'
#' ge_list to flat file. Used when samples are in different files on GEO
#'
#' @return path to raw, prepped input files
#'
#' @param ge_list list of gene expression tables (matrix, data.frame, etc)
#' @param supp_files_dir path to base directory (via \code{.get_supp_files_dir()})
#' @param study study accession eg \code{SDY269}
.ge_list_to_flat_file <- function(ge_list, supp_files_dir, study){
  ge_df <- Reduce(f = function(x, y) {merge(x, y)}, ge_list)
  input_files <- .write_raw_expression(ge_df, supp_files_dir, study)
}
#######################################
###            MAPPING              ###
#######################################

#' Standardize probe column name
#'
#' Different platforms use different names for the feature id column.
#' This will update this column to be called \code{feature_id}
#'
#' @param exprs_dt data.table of gene expression with one column per sample,
#' one row per feature.
#'
.map_feature_id_col <- function(exprs_dt){
  if ( !any(grepl("feature_id", colnames(exprs_dt))) ) {

    # If Illumina from Immport
    prbCol <- grep("id_ref", colnames(exprs_dt), ignore.case = TRUE)

    # If RNAseq then accept gene* or V1 col
    if (length(prbCol) == 0) {
      prbCol <- grep("gene|^V1$", colnames(exprs_dt), ignore.case = TRUE)
    }

    # In case of features in rownames, e.g. from GEO
    if (length(prbCol) == 0) {
      prbCol <- "rn"
    }

    setnames(exprs_dt, prbCol, "feature_id")
  }
  return(exprs_dt)
}


#' Summarize Matrix
#'
#' Summarize any duplicated genes so that there is
#' one entry per gene.
#'
#' @details This function summarizes normalized expression values by mean
#' into a table with one entry per gene name based on \code{feature_gene_map}.
#' The input data should be in log2 space.
#'
#' @param ge_dt data.table of gene expression values. Features should be listed in
#' \code{feature_id} column, with additional columns for each sample.
#' @param feature_gene_map data.table mapping original gene id (eg probe id, ensembl id, gene alias)
#' to gene symbol. It should have two columns: \code{featureid} (original) and \code{genesymbol} (updated).
#'
#' @export
summarize_by_gene_symbol <- function(ge_dt,
                             feature_gene_map){
  em <- copy(mx)
  em[ , gene_symbol := feature_gene_map[match(em$feature_id, feature_gene_map$featureid), genesymbol] ]
  em <- em[ !is.na(gene_symbol) & gene_symbol != "NA" ]
  sum_exprs <- em[ , lapply(.SD, mean),
                   by = "gene_symbol",
                   .SDcols = grep("^BS", colnames(em)) ]
}

#' Write matrix to file system
#'
#' @details Writes four versions of flat tsv files:
#'   1. \code{<matrix_name>.raw.tsv}: raw, background-corrected values
#'   2. \code{<matrix_name>.tsv}: normalized values
#'   3. \code{<matrix_name>.summary.tsv}: normalized values, summarized by gene symbol
#'   (based on current annotation)
#'   4. \code{<matrix_name>.summary.orig}: normalized values, summarized by gene symbol
#'   (based on original annotation)
#'
#' @param pipeline.root pipeline.root
#' @param matrix_name Name of matrix
#' @param exprs data.table of raw (background-corrected) expression (from \code{makeRawMatrix()})
#' @param norm_exprs data.table of normalized expression (from \code{normalizeMatrix()})
#' @param sum_exprs data.table of summarized expression (from \code{summarizeMatrix()})
#' @param verbose Write verbose print statements?
#'
#' @export
write_matrix <- function(output_dir,
                         matrix_name,
                         exprs,
                         norm_exprs,
                         sum_exprs,
                         verbose = FALSE){

  .write_em <- function(df, file_name, verbose){
    if ( verbose ) message("Writing ", file_name, "...")
    fwrite(df,
           file = file.path(output_dir,
                            file_name),
           sep = "\t",
           quote = FALSE,
           row.names = FALSE)
  }

  # Raw EM
  .write_em(exprs, paste0(matrix_name, ".tsv.raw"), verbose)

  # Normalized EM
  .write_em(norm_exprs, paste0(matrix_name, ".tsv"), verbose)

  # summary EM
  .write_em(sum_exprs, paste0(matrix_name, ".tsv.summary"), verbose)

  # original summary EM assuming run created with _orig FasId
  .write_em(sum_exprs, paste0(matrix_name, ".tsv.summary.orig"), verbose)

}


#' Make Raw Matrix
#'
#' @return A data.table containing raw counts or background-corrected probe
#' intensities with a feature_id column and one column per biosample_accession
#'
#' @param meta_data list of study-specific meta data
#' @param study study accession eg \code{SDY269}
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
#' @param input_files input file names
#' @export
make_raw_matrix <- function(meta_data,
                            gef,
                            study,
                            input_files,
                            analysis_dir) {

  # Get raw files from GEO or prepare ImmPort flat files
  # At end of this step, there should be a single "cohort_type_raw_expression.txt"
  # file for non-affymetrix studies
  if ( meta_data$file_location %in% c("gsm_soft", "gsm_supp_files", "gse_supp_files") ) {
    input_files <- .prep_geo_files(study,
                                   gef,
                                   meta_data,
                                   input_files,
                                   analysis_dir)
  } else {
    input_files <- .prep_immport_files(study,
                                       gef,
                                       platform,
                                       input_files,
                                       analysis_dir)
  }

  # Generate background corrected raw matrices for affy and illumina
  # For RNAseq pass through raw counts file.
  exprs_dt <- switch(
    meta_data$platform,
    "Affymetrix" = .process_affy(input_files),
    "Illumina" = .process_illumina(input_files),
    "NA" = .process_rna_seq(input_files),
    fread(input_files)
  )



  # Ensure probe column is named 'feature_id'
  exprs_dt <- .map_feature_id_col(exprs_dt)

  # Ensure all probes have names. Not a problem for most studies.
  # Note that values can still be NA here and may be due to a handful
  # of samples having problems (e.g. SDY224 / SDY212).  These NA values
  # are removed during normalization.
  exprs_dt <- exprs_dt[ !is.na(feature_id) & feature_id != "" ]

  # Map expsample or geo accession to biosample accession
  exprs_dt <- .map_sampleid_to_biosample_accession(exprs_dt, gef)

  # Subset to biosamples in selected_biosamples from UI
  exprs_dt <- exprs_dt[ , colnames(exprs_dt) %in% c("feature_id", gef$biosample_accession),
                        with = FALSE ]

  # Check that all gef samples are in the matrix - may not be the case if some
  # were removed for QC issues.  User should re-run matrix with these samples
  # removed.
  if ( !all(gef$biosample_accession %in% colnames(exprs_dt)) ) {
    stop("Some selected biosamples are not found in raw-matrix.",
         "These may have been removed for QC issues. Please check and re-run")
  }

  # Ensure colOrder
  sampleids <- grep("feature_id", names(exprs_dt), invert = TRUE)
  feature_id <- grep("feature_id", names(exprs_dt))
  setcolorder(exprs_dt, c(feature_id, sampleids))

  return(exprs_dt)
}


