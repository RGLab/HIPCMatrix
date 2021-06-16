
#-------------------------------
# PIPELINE WRAPPER
#-------------------------------

# onCL means onCommandLine and avoids writing out extra

#' Run Create Matrix
#'
#' Runs all steps of creating matrix in ImmuneSpace.
#'
#' @param labkey.url.base labkey.url.base
#' @param study Study accession
#' @param pipeline.root pipeline.root
#' @param analysis.directory analysis.directory
#' @param output.tsv Name of output matrix
#' @param selectedBiosamples List of biosamples (should be from one
#' cell type from one cohort)
#' @param fasId Feature Annotation Set ID
#' @param taskOutputParams taskOutputParams
#' @param onCL running on command line? (T/F)
#'
#' @import ImmuneSpaceR
#' @import Rlabkey
#' @import data.table
#'
#' @export
runCreateMx <- function(labkey.url.base,
                        study,
                        base_dir,
                        analysis.directory,
                        output.tsv,
                        selectedBiosamples,
                        fasId,
                        taskOutputParams,
                        onCL = FALSE) {

  # -------------------------------- RETRIEVE INPUTS ----------------------------------
  # For printing and con
  stopifnot(grepl("^SDY\\d+$", study))
  matrix_name <- gsub(".tsv", "", output.tsv)

  if (onCL == TRUE) {
    print(paste(study, matrix_name))
  }

  # Check that output filepath exists before starting run
  outPath <- file.path(pipeline.root, "analysis/exprs_matrices")
  if (!dir.exists(outPath)) {
    dir.create(outPath)
  }

  # Check that feature2gene mapping is available prior to doing work
  co <- labkey.setCurlOptions(ssl_verifyhost = 2, sslversion = 1)

  FAS_filter <- makeFilter(c("FeatureAnnotationSetId/RowId",
                             "IN",
                             fasId))

  feature_gene_map <- data.table(labkey.selectRows(baseUrl = labkey.url.base,
                                                   folderPath = "/Studies/",
                                                   schemaName = "Microarray",
                                                   queryName = "FeatureAnnotation",
                                                   colFilter = FAS_filter,
                                                   colNameOpt = "rname",
                                                   colSelect = c("featureid","genesymbol")))

  if (nrow(feature_gene_map) == 0) {
    stop("The downloaded feature annotation set has 0 rows.")
  }

  # Specifying study and onTest so that code can be used in either module / pipeline,
  # which currently pulls lub and lup from javascript calls, or CL work.
  onTest <- labkey.url.base == "https://test.immunespace.org"
  con <- CreateConnection(study = study, onTest = onTest)

  # Create GEF
  bs_filter <- makeFilter(c("biosample_accession", "IN", gsub(",", ";", selectedBiosamples)))
  gef <- con$getDataset("gene_expression_files",
                        colFilter = bs_filter,
                        original_view = TRUE,
                        reload = TRUE)

  # ensure single cohort for processing
  if (length(unique(gef$arm_name)) > 1) {
    stop("There are more than one cohort selected in this HIPCMatrix run.")
  }

  # ensure each expsample has unique biosample, otherwise summarization is thrown off
  if (length(gef$biosample_accession) != length(gef$expsample_accession)) {
    stop("Experiment samples do not have unique biosample accessions.")
  }

  metaData <- getMetaData(study = study,
                          gef = gef,
                          fasId = fasId)

  # ----------------------------- PROCESSING -------------------------------------

  # Identify correct input_files
  if ( metaData$file_location == "immport" ) {

    input_files <- gef$file_info_name[ grep("cel$|txt$|tsv$|csv$",
                                            gef$file_info_name,
                                            ignore.case = TRUE)]
    input_files <- paste0(pipeline.root, "/rawdata/gene_expression/", input_files)
    input_files <- unique(input_files[ file.exists(input_files) ])

  } else if ( metaData$file_location == "custom" ) {

    suffix <- ifelse(metaData$specialCase, "author_data/", "raw_counts/")
    sdyGEDir <- paste0("/share/files/Studies/", study, "/@files/rawdata/gene_expression/")
    rawDir <- paste0(sdyGEDir, suffix)

    input_files <- list.files(rawDir, full.names = TRUE)

    uniqueSubString <- ifelse(specialCase, "GA", "Header")
    input_files <- input_files[ grep("uniqueSubString", input_files, invert = TRUE) ]

  } else {

    input_files <- NA
    gef <- gef[ !is.na(gef$geo_accession), ]

  }

  # Create three versions of matrix
  exprs <- make_raw_matrix(metaData = metaData,
                           gef = gef,
                           study = study,
                           input_files = input_files)

  norm_exprs <- normalize_matrix(exprs, metaData)

  sum_exprs <- summarize_matrix(norm_exprs, feature_gene_map)

  # ------------------------------ OUTPUT ------------------------------------------
  writeMatrix(pipeline.root, output.tsv, exprs, norm_exprs, sum_exprs, onCL)

  # This file gets cleaned up anyway, so not worrying about it onCL
  if (onCL == FALSE) {
    outProps = file(description = taskOutputParams, open = "w")
    cat(file = outProps, sep="", "name\tvalue\n")
    cat(file = outProps, sep="", "assay run property, cohort\t", unique(gef$cohort), "\n")
    flush(con = outProps)
    close(con = outProps)
  }

  # create copy of CM.R script from run time, after checking to be sure analysis
  # directory is in place. It is missing from some studies for some reason.
  if (!dir.exists(analysis.directory)) {
    dir.create(analysis.directory, recursive = TRUE)
  }

  # Allow for work on server or local
  LKModules <- "~/LabKeyModules"

  file.copy(from = file.path(LKModules, "HIPCMatrix/pipeline/tasks/create-matrix.R"),
            to = paste0(analysis.directory, "/", output.tsv, "-create-matrix-snapshot.R"))

  file.copy(from = file.path(LKModules, "HIPCMatrix/pipeline/tasks/runCreateMx.R"),
            to = paste0(analysis.directory, "/", output.tsv, "-runCM-snapshot.R"))

  # write out tsv of vars to make later replication of results easier
  varDf <- data.frame(labkey.url.base = labkey.url.base,
                      labkey.url.path = labkey.url.path,
                      pipeline.root = pipeline.root,
                      analysis.directory = analysis.directory,
                      selectedBiosamples = selectedBiosamples,
                      fasId = fasId,
                      taskOutputParams = taskOutputParams,
                      output.tsv = output.tsv,
                      stringsAsFactors = FALSE)


  write.table(varDf,
              file = paste0(analysis.directory, "/create-matrix-vars.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}



#' Make Raw Matrix
#'
#' @return A data.table containing raw counts or background-corrected probe
#' intensities with a feature_id column and one column per biosample_accession
#'
#' @param metaData list of study-specific meta data
#' @param study study accession eg \code{SDY269}
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
#' @param input_files input file names
#' @export
make_raw_matrix <- function(metaData,
                            gef,
                            study,
                            input_files) {

  # Get raw files from GEO or prepare ImmPort flat files
  # At end of this step, there should be a single "cohort_type_raw_expression.txt"
  # file for non-affymetrix studies
  if ( metaData$file_location %in% c("gsm_soft", "gsm_supp_files", "gse_supp_files") ) {
    input_files <- .prep_geo_files(study, gef, metaData, input_files)
  } else {
    input_files <- .prep_immport_files(study, gef, metaData, input_files)
  }

  # Generate background corrected raw matrices for affy and illumina
  # For RNAseq pass through raw counts file.
  if ( metaData$platform == "Affymetrix" ) {
    exprs_dt <- .process_affy(input_files, gef, metaData)
  } else {
    if ( metaData$platform == "Illumina" ) {
      exprs_dt <- .process_illumina(input_files)
    } else if (metaData$platform == "NA"){
      exprs_dt <- .process_rna_seq(input_files, study)
    } else {
      exprs_dt <- fread(input_files)
    }

    # Ensure probe col is named 'feature_id'
    exprs_dt <- .map_feature_id_col(exprs_dt)

    # Ensure all probes have names. Not a problem for most studies.
    # Note that values can still be NA here and may be due to a handful
    # of samples having problems (e.g. SDY224 / SDY212).  These NA values
    # are removed during normalization.
    exprs_dt <- exprs_dt[ !is.na(feature_id) & feature_id != "" ]
  }

  # Map expsample or geo accession to biosample accession
  exprs_dt <- .map_sampleid_to_biosample_accession(exprs_dt, gef)

  # Subset to biosamples in selectedBiosamples from UI
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
  setcolorder(exprs_dt, c(featur_id, sampleids))

  return(exprs_dt)
}


#' Write matrix to file system
#'
#' @details Writes four verisons of flat tsv files:
#'   1. \code{<matrix_name>.raw.tsv}: raw, background-corrected values
#'   2. \code{<matrix_name>.tsv}: normalized values
#'   3. \code{<matrix_name>.summary.tsv}: normalized values, summarized by gene symbol
#'   (based on current annotation)
#'   4. \code{<matrix_name>.summary.orig}: normalized values, summarized by gene symbol
#'   (based on original annotation)
#'
#' @param pipeline.root pipeline.root
#' @param output.tsv Name of matrix
#' @param exprs data.table of raw (background-corrected) expression (from \code{makeRawMatrix()})
#' @param norm_exprs data.table of normalized expression (from \code{normalizeMatrix()})
#' @param sum_exprs data.table of summarized expression (from \code{summarizeMatrix()})
#' @param onCL running on command line? (T/F)
#'
#' @export
write_matrix <- function(pipeline.root,
                         output.tsv,
                         exprs,
                         norm_exprs,
                         sum_exprs,
                         onCL = FALSE){

  .writeTbl <- function(df, onCL, baseNm){
    dir <- ifelse(onCL == TRUE, paste0(pipeline.root, "/analysis/exprs_matrices/"), "")
    write.table(df,
                file = paste0(dir, baseNm),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
  }

  # Raw EM
  .writeTbl(exprs, onCL, paste0(output.tsv, ".raw"))

  # Normalized EM
  .writeTbl(norm_exprs, onCL, output.tsv)

  # summary EM
  .writeTbl(sum_exprs, onCL, paste0(output.tsv, ".summary"))

  # original summary EM assuming run created with _orig FasId
  .writeTbl(sum_exprs, onCL, paste0(output.tsv, ".summary.orig"))

}
