#' Get Base Directory
#'
#' Create run-specific base directory for supplementary files.
#' Important that run specific directory is created to avoid overlap with previous runs.
#' This done with unique cell_type * arm_accession string (aka cohort_type)
#'
#' @param study study accession eg \code{SDY269}
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
.get_supp_files_dir <- function(study, gef){
  supp_files_dir <- file.path("/share/files/Studies",
                    study,
                    "@files/rawdata/gene_expression/supp_files",
                    paste0(unique(gef$type),
                           "_",
                           unique(gef$arm_accession)))
  dir.create(supp_files_dir, recursive = TRUE)
  return(supp_files_dir)
}

#-------------------------------
# CREATE MATRIX FN
#-------------------------------

#' Make Raw Matrix
#'
#' @return A data.table with a feature_id column and one column per biosample_accession
#'
#' @param metaData list of study-specific meta data
#' @param study study accession eg \code{SDY269}
#' @param gef result of ISCon$getDataset("gene_expression_files") for one run.
#' @param input_files input file names
#' @export
makeRawMatrix <- function(metaData,
                          gef,
                          study,
                          input_files) {

  # Get raw files from GEO or prepare ImmPort flat files
  # At end of this step, there should be a single "cohort_type_raw_expression.txt"
  # file for non-affymetrix studies
  if ( metaData$file_location %in% c("gsm_soft", "gsm_supp_files", "gse_supp_files") ) {
    input_files <- .prep_geo_files(study, gef, metaData, input_files)
  } else {
    input_files <- .prepImmportFls(study, gef, metaData, input_files)
  }

  # Generate background corrected raw matrices for affy and illumina
  # For RNAseq pass through raw counts file.
  if ( metaData$platform == "Affymetrix" ) {
    exprs <- .processAffy(input_files, gef, metaData)
  } else {
    if ( metaData$platform == "Illumina" ) {
      exprs <- .processIllumina(input_files)
    } else if (metaData$platform == "NA"){
      exprs <- .processRnaSeq(input_files, study)
    } else {
      exprs <- fread(input_files)
    }

    # Ensure probe col is named 'feature_id'
    exprs <- .mapFeatureIdCol(exprs)

    # Ensure all probes have names. Not a problem for most studies.
    # Note that values can still be NA here and may be due to a handful
    # of samples having problems (e.g. SDY224 / SDY212).  These NA values
    # are removed during normalization.
    exprs <- exprs[ !is.na(feature_id) & feature_id != "" ]
  }

  # Map expsample or geo accession to biosample accession
  exprs <- .mapAccToBs(exprs, gef)

  # Subset to biosamples in selectedBiosamples from UI
  exprs <- exprs[ , colnames(exprs) %in% c("feature_id", gef$biosample_accession),
                  with = FALSE ]

  # Check that all gef samples are in the matrix - may not be the case if some
  # were removed for QC issues.  User should re-run matrix with these samples
  # removed.
  if (!all(gef$biosample_accession %in% colnames(exprs))) {
    stop("Some selected biosamples are not found in raw-matrix. These may have been removed for QC issues. Please check and re-run")
  }

  # Ensure colOrder
  smpls <- grep("feature_id", names(exprs), invert = TRUE)
  fid <- grep("feature_id", names(exprs))
  setcolorder(exprs, c(fid, smpls))

  return(exprs)
}

#' Normalize Matrix
#'
#'
#' @param exprs data.table of expression
#' @param study study accession eg \code{SDY269}
#' @param metaData list of study-specific meta data
#'
#' @import DESeq2
#'
#' @export
normalizeMatrix <- function(exprs, study, metaData){
  # data.table so MUST COPY to prevent changes in later work
  em <- copy(exprs)

  # Oct 2018 - following studies have some NA values.
  # SDY212 - Known issue and documented in ImmPort Jira
  # SDY1289 - Due to single cohort having multiple batches with different anno.
  # SDY224 - Controls
  em <- em[ complete.cases(em), ]

  # already processed and raw FASTQ data only in SRA, no raw count matrix easily available
  if (metaData$noRaw == TRUE) {
    return(em)
  }

  rnames <- em[ , feature_id ]
  em[ , feature_id := NULL ]

  # Must ensure numeric values as conversion back to character can happen with
  # casting as matrix.
  em <- as.matrix(apply(em, 2, as.numeric))

  # no platform  == isRNASeq
  if (metaData$platform == "NA") {
    # newCountDataSet does not take duplicated column names, so assign temporary unique names
    orginal_colnames <- colnames(em)
    colnames(em) <- seq_len(ncol(em))

    cds <- DESeq2::DESeqDataSetFromMatrix(countData = em,
                                          colData = data.frame(cohort = rep(1, ncol(em))),
                                          design = ~ 1)
    cds <- DESeq2::estimateSizeFactors(cds)
    cdsBlind <- DESeq2::estimateDispersions(cds)
    vsd <- DESeq2::varianceStabilizingTransformation(cdsBlind)
    norm_exprs <- SummarizedExperiment::assay(vsd)
    colnames(norm_exprs) <- orginal_colnames
  } else {
    cnames <- colnames(em)
    norm_exprs <- preprocessCore::normalize.quantiles(em)
    colnames(norm_exprs) <- cnames
    norm_exprs <- pmax(norm_exprs, 1)
    if (max(norm_exprs) > 100) {
      norm_exprs <- log2(norm_exprs)
    }
  }

  norm_exprs <- data.table(norm_exprs)
  norm_exprs[ , feature_id := rnames ]

  smpls <- grep("feature_id", names(norm_exprs), invert = TRUE)
  fid <- grep("feature_id", names(norm_exprs))
  setcolorder(norm_exprs, c(fid, smpls))

  return(norm_exprs)
}

#' Summarize Matrix
#'
#' Summarize any duplicated genes so that there is
#' one entry per gene.
#'
#' @details This function summarizes normalized expression values into a table
#' with one entry per gene name by:
#'   1. Map original feature names to gene symbols based on \code{f2g} table
#'   2. Remove any features that did not map to a gene symbol
#'   3. Summarize duplicated gene symbols by mean of the (normalized) values.
#'
#' @param norm_exprs data.table of normalized gene expression values
#' @param f2g data.table mapping original gene id (eg probe id, ensembl id, gene alias)
#' to gene symbol. It should have two columns: \code{featureid} (original) and \code{genesymbol} (updated).
#'
#' @export
summarizeMatrix <- function(norm_exprs, f2g){
  em <- copy(norm_exprs)
  em[ , gene_symbol := f2g[match(em$feature_id, f2g$featureid), genesymbol] ]
  em <- em[ !is.na(gene_symbol) & gene_symbol != "NA" ]
  sum_exprs <- em[ , lapply(.SD, mean),
                   by = "gene_symbol",
                   .SDcols = grep("^BS", colnames(em)) ]
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
writeMatrix <- function(pipeline.root,
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
                        pipeline.root,
                        analysis.directory,
                        output.tsv,
                        selectedBiosamples,
                        fasId,
                        taskOutputParams,
                        onCL = FALSE){

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

  f2g <- data.table(labkey.selectRows(baseUrl = labkey.url.base,
                                      folderPath = "/Studies/",
                                      schemaName = "Microarray",
                                      queryName = "FeatureAnnotation",
                                      colFilter = FAS_filter,
                                      colNameOpt = "rname",
                                      colSelect = c("featureid","genesymbol")))

  if (nrow(f2g) == 0) {
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
  exprs <- makeRawMatrix(metaData = metaData,
                         gef = gef,
                         study = study,
                         input_files = input_files)

  norm_exprs <- normalizeMatrix(exprs, study, metaData)

  sum_exprs <- summarizeMatrix(norm_exprs, f2g)

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
                      stringsAsFactors = FALSE
  )

  write.table(varDf,
              file = paste0(analysis.directory, "/create-matrix-vars.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}
