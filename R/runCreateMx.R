
#-------------------------------
# Run EVERYTHING
#-------------------------------

# onCL means onCommandLine and avoids writing out extra

#' Run Create Matrix
#'
#' Runs all steps of creating matrix in ImmuneSpace.
#'
#' @param labkey.url.base labkey.url.base
#' @param study Study accession
#' @param matrix_name Name of the matrix
#' @param base_dir Base directory for all analysis and output files.
#'  If \code{output_dir}, \code{analysis_dir}, \code{debug_dir} are not
#'  defined, they will be created under \code{base_dir}. When running in
#'  the HIPCMatrix module on ImmuneSpace servers,
#'  this should be the same as \code{pipeline.root}.
#' @param output_dir Directory to write final matrices. On ImmuneSpace servers:
#' `/share/files/Studies/<study_accession>/@files/analysis/exprs_matrices`
#' @param analysis_dir Directory to write intermediate files. On servers,
#'  `/share/files/Studies/<study_accession>/@files/rawdata/gene_expression`
#' @param debug_dir Directory to write debug files. On servers,
#' it will be the same as \code{analysis.directory} in the HIPCMatrix pipeline,
#' @param selected_biosamples String listing biosample accessions to include,
#' separated by commas (should be from one cell type from one cohort)
#' @param fas_id Feature Annotation Set ID
#' @param taskOutputParams taskOutputParams
#' @param verbose Print verbose debug statements?
#' @param snapshot write copies of source files to `debug_dir`?
#' @param reload Force re-download of supplementary files from GEO?
#'
#' @import ImmuneSpaceR
#' @import Rlabkey
#' @import data.table
#'
#' @export
runCreateMx <- function(study,
                        matrix_name,
                        selected_biosamples,
                        fas_id,
                        labkey.url.base = "https://www.immunespace.org/",
                        base_dir = file.path(
                          "/share",
                          "files",
                          "Studies",
                          study,
                          "@files"
                        ),
                        output_dir = file.path(
                          base_dir,
                          "analysis",
                          "exprs_matrices"
                        ),
                        analysis_dir = file.path(
                          base_dir,
                          "rawdata",
                          "gene_expression"
                        ),
                        debug_dir = file.path(
                          analysis_dir,
                          "create-matrix",
                          matrix_name
                        ),
                        taskOutputParams = NULL,
                        verbose = TRUE,
                        snapshot = TRUE,
                        reload = TRUE) {
  if (verbose) {
    log_message(
      "Running runCreateMx using HIPCMatrix version ",
      utils::packageVersion("HIPCMatrix")
    )
  }
  # -------------------------------- RETRIEVE INPUTS ----------------------------------
  # For printing and con
  stopifnot(grepl("^SDY\\d+$", study))

  if (verbose) {
    log_message(matrix_name)
  }

  # Check that output filepath exists before starting run
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Check that feature2gene mapping is available prior to doing work
  # TODO: Why?
  co <- labkey.setCurlOptions(ssl_verifyhost = 2, sslversion = 1)

  FAS_filter <- makeFilter(c(
    "FeatureAnnotationSetId/RowId",
    "IN",
    fas_id
  ))

  feature_gene_map <- data.table(labkey.selectRows(
    baseUrl = labkey.url.base,
    folderPath = "/Studies/",
    schemaName = "Microarray",
    queryName = "FeatureAnnotation",
    colFilter = FAS_filter,
    colNameOpt = "rname",
    colSelect = c("featureid", "genesymbol")
  ))

  if (nrow(feature_gene_map) == 0) {
    stop("The downloaded feature annotation set has 0 rows.")
  }

  # Specifying study and onTest so that code can be used in either module / pipeline,
  # which currently pulls lub and lup from javascript calls, or CL work.
  onTest <- labkey.url.base == "https://test.immunespace.org"
  con <- CreateConnection(study = study, onTest = onTest)

  # Create GEF
  bs_filter <- makeFilter(c("biosample_accession", "IN", gsub(",", ";", selected_biosamples)))
  gef <- con$getDataset("gene_expression_files",
    colFilter = bs_filter,
    original_view = TRUE,
    reload = TRUE
  )

  # ensure single cohort for processing
  if (length(unique(gef$arm_name)) > 1) {
    stop("There are more than one cohort selected in this HIPCMatrix run.")
  }

  # ensure each expsample has unique biosample, otherwise summarization is thrown off
  if (length(gef$biosample_accession) != length(gef$expsample_accession)) {
    stop("Experiment samples do not have unique biosample accessions.")
  }

  meta_data <- get_meta_data(
    study = study,
    gef = gef,
    fas_id = fas_id,
    baseUrl = con$config$labkey.url.base
  )

  # ----------------------------- PROCESSING -------------------------------------


  input_files <- retrieve_input_files(
    study = study,
    gef = gef,
    meta_data = meta_data,
    analysis_dir = analysis_dir,
    verbose = verbose,
    reload = reload
  )

  # Create three versions of matrix
  exprs <- make_raw_matrix(
    platform = meta_data$platform,
    gef = gef,
    input_files = input_files,
    verbose = verbose
  )

  norm_exprs <- normalize_matrix(
    exprs,
    meta_data$platform,
    verbose = verbose
  )

  sum_exprs <- summarize_by_gene_symbol(
    norm_exprs,
    feature_gene_map,
    verbose = verbose
  )

  # ------------------------------ OUTPUT ------------------------------------------
  write_matrix(
    output_dir = output_dir,
    matrix_name = matrix_name,
    exprs = exprs,
    norm_exprs = norm_exprs,
    sum_exprs = sum_exprs,
    verbose = verbose
  )

  if (!is.null(taskOutputParams)) {
    outProps <- file(description = taskOutputParams, open = "w")
    cat(file = outProps, sep = "", "name\tvalue\n")
    cat(file = outProps, sep = "", "assay run property, cohort\t", unique(gef$cohort), "\n")
    cat(file = outProps, sep = "", "assay run property, version\t", as.character(packageVersion("HIPCMatrix")), "\n" )
    cat(file = outProps, sep = "", "assay run property, hash\t", as.character(sessioninfo:::pkg_desc("HIPCMatrix")$GithubSHA1), "\n" )
    flush(con = outProps)
    close(con = outProps)
  }

  if (snapshot) {

    if (verbose) log_message("Writing snapshot of sources to ", debug_dir)
    # create copy of CM.R script from run time, after checking to be sure analysis
    # directory is in place. It is missing from some studies for some reason.
    if (!dir.exists(debug_dir)) {
      dir.create(debug_dir, recursive = TRUE)
    }

    # Allow for work on server or local
    LKModules <- "/labkey/git/LabKeyModules"

    file.copy(
      from = file.path(LKModules, "HIPCMatrix/pipeline/tasks/create-matrix.R"),
      to = file.path(debug_dir, paste0(matrix_name, "-create-matrix-snapshot.R"))
    )

    writeLines(
      as.character(utils::packageVersion("HIPCMatrix")),
      file.path(debug_dir, paste0(matrix_name, "-HIPCMatrix-version.txt"))
    )

    # write out tsv of vars to make later replication of results easier
    varDf <- data.table(
      study = study,
      matrix_name = matrix_name,
      selected_biosamples = selected_biosamples,
      fas_id = fas_id,
      labkey.url.base = labkey.url.base,
      base_dir = base_dir,
      output_dir = output_dir,
      analysis_dir = analysis_dir,
      debug_dir = debug_dir,
      taskOutputParams = taskOutputParams,
      verbose = verbose,
      snapshot = snapshot
    )


    fwrite(varDf,
      file = file.path(debug_dir, "create-matrix-vars.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}
