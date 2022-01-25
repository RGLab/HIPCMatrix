#' @include HMX.R
NULL

#' Find Differentially Expressed genes using Empirical Bayes Statistics
#'
#' Compute differential expression over time on microarray data using
#' empirical Bayes statistics on a linear model using the
#' \code{\link[limma]{eBayes}} method from the \code{limma} package. All
#' differentially expressed genes for each timepoint are reported.
#'
#' @return list with one entry for each non-baseline timepoint
#'
#' @param eset expressionset object with normalized microarray expression.
#' @param p_val_cutoff adjusted p-value cutoff for determining differentially
#' expressed genes.
#'
#' @export
find_de_genes_eBayes <- function(eset,
                                 p_val_cutoff = 0.02) {
  contrast <- c("study_time_collected", "study_time_collected_unit")
  pd <- data.table(Biobase::pData(eset))
  pd <- pd[, coef := do.call(paste, .SD), .SDcols = contrast]
  if (length(unique(pd$coef)) < 2) {
    stop("<2 groups found. Cannot perform DE analysis. ")
  }

  to_drop <- unique(pd[study_time_collected <= 0, coef])
  pd <- pd[coef %in% to_drop, coef := "baseline"]
  tmp <- grep("baseline", value = TRUE, invert = TRUE, gtools::mixedsort(unique(pd$coef)))
  pd <- pd[, coef := factor(coef, levels = c("baseline", tmp))] # preps coef col for use in model

  if (length(unique(pd$participant_id)) < 2) {
    stop("Fewer than 2 participants with multiple timepoints! ")
  }

  mm <- stats::model.matrix(stats::formula("~participant_id + coef"), pd)

  if (dim(mm)[[1]] < dim(mm)[[2]]) {
    stop("Not enough subjects to perform analysis")
  }
  # NOTE:  ES is the normalized expressionset

  # TODO:  add more notes here and on notion
  fit <- limma::lmFit(eset, mm)
  fit <- tryCatch(
    limma::eBayes(fit),
    error = function(e) {
      return(e)
    }
  )

  if (!is.null(fit$message)) {
    stop("Linear model not able to be fit: ", fit$message)
  }

  coefs <- grep("^coef", colnames(mm), value = TRUE)
  tt_list <- list()
  for (coef in coefs) {

    # Check that coef can be used
    tt <- data.table(limma::topTable(fit, coef = coef, number = Inf))
    if (all(is.na(tt$adj.P.Val))) {
      log_message(coef, " has all NA values for adj.P.Val. Skipping to next coef.")
      next()
    }

    tt <- if (sum(tt$adj.P.Val < p_val_cutoff, na.rm = T) < 100) {
      tt[order(adj.P.Val)][1:min(nrow(tt), 100)]
    } else {
      tt[adj.P.Val < p_val_cutoff]
    }

    tt[, coefficient := gsub("coef", "", coef)]
    tt_list[[coef]] <- tt
  }

  tt_list
}


##### ----------------- HMX methods ------------------------ #####

# See HMX$run_de_analysis
run_de_analysis <- function(con, rerun = FALSE) {
  if (grepl("^IS", con$study)) {
    stop("run_de_analysis is not designed for ImmuneSignatures!")
  }

  if (!grepl("SDY", con$study)) {
    stop("run run_de_analysis one study at a time!")
  }

  if ("de_results" %in% names(con$cache) & !rerun) {
    return(con$cache$de_results)
  }

  if (con$config$verbose) {
    log_message("Running Differential Expression Analysis...")
  }

  de_runs <- con$get_de_compatible_runs()

  if (nrow(de_runs) == 0) {
    log_message("Insufficient data for all cohorts. Analysis not run.")
    con$cache[["de_runs"]] <- NULL
    con$cache[["de_results"]] <- NULL
    return(NULL)
  }
  if (con$config$verbose) {
    log_message(nrow(de_runs), " cohorts/timepoints found for differential expression")
  }

  GEA_list <- vector("list")
  GEAR_list <- vector("list")

  runs <- con$cache$GE_matrices$name
  gef <- con$getDataset("gene_expression_files", original_view = T)

  idx <- 1 # analysis accession key
  for (run in runs) {
    eset <- con$getGEMatrix(matrixName = run, outputType = "normalized", annotation = "latest")
    tt_list <- tryCatch(find_de_genes_eBayes(eset),
      error = function(e) e
    )

    if ("error" %in% class(tt_list)) {
      log_message("Could not run DE analysis for ", run, ": ", tt_list$message)
      next()
    }

    cm <- con$getDataset("cohort_membership")
    cm <- unique(cm[, list(cohort, arm_accession)])

    # subset GEF to match run samples
    gefSub <- gef[gef$biosample_accession %in% eset$biosample_accession, ]

    # create cohort_type identifier so unique to runs, which are cohort * cell_type
    cohort_type <- unique(paste(gefSub$cohort, gefSub$type, sep = "_"))

    arm_accession <- cm[cohort == unique(eset$cohort), arm_accession]

    for (coef in names(tt_list)) {
      # One row in GEA per coef per run

      analysis_accession <- paste0("GEA", idx)
      # arm_name[ is.null(arm_name) ] <- NA
      description <- paste0("Differential expression in ", run, ", ", gsub("^coef", "", coef), " vs. baseline")


      GEA_list[[idx]] <- data.table(
        analysis_accession = analysis_accession,
        expression_matrix = run,
        arm_name = cohort_type,
        arm_accession = arm_accession,
        coefficient = gsub("^coef", "", coef),
        description = description
      )

      tt <- tt_list[[coef]]
      if (nrow(tt) > 0) {
        tt[, analysis_accession := analysis_accession]
        GEAR_list[[idx]] <- data.table(tt)
      }

      idx <- idx + 1
    }
  }

  de_runs <- rbindlist(GEA_list)
  de_results <- rbindlist(GEAR_list)

  setnames(
    de_results,
    c("FeatureId", "gene_symbol", "adj.P.Val", "AveExpr", "logFC", "P.Value", "t"),
    c("feature_id", "gene_symbol", "adj_p_val", "ave_expr", "log_fc", "p_value", "statistic")
  )

  setcolorder(
    de_results,
    neworder = c(
      "feature_id",
      "gene_symbol",
      "log_fc",
      "ave_expr",
      "statistic",
      "p_value",
      "adj_p_val",
      "B",
      "analysis_accession",
      "coefficient"
    )
  )

  con$cache[["de_runs"]] <- de_runs
  con$cache[["de_results"]] <- de_results
  return(de_results)
}

# see HMX$upload_de_analysis_results
upload_de_analysis_results <- function(con) {
  if (!grepl("SDY", con$study)) {
    stop("Please run for one study at a time.")
  }
  de_results <- con$run_de_analysis()

  de_runs <- con$cache$de_runs
  if (is.null(de_results)) {
    log_message(
      "Differential expression results not found. ",
      "Skipping upload."
    )
    return(invisible(con))
  }

  if (nrow(de_results) == 0) {
    stop("No differentially expressed genes found.")
  }

  if (nrow(de_runs) > 0) {

    # delete old GEA
    currGEA <- tryCatch(
      ImmuneSpaceR:::.getLKtbl(
        con = con,
        schema = "gene_expression",
        query = "gene_expression_analysis"
      ),
      error = function(e) {
        return(e)
      }
    )

    if (!is.null(currGEA$message)) {
      log_message(paste0("Error: ", currGEA$message))
      stop("May need to turn on DEA module to allow gene_expression schema.")
    }

    log_message("Deleting all ", nrow(currGEA), " rows in gene_expression_analysis table...")
    if (nrow(currGEA) != 0) {
      deleteGEA <- labkey.deleteRows(
        baseUrl = con$config$labkey.url.base,
        folderPath = con$config$labkey.url.path,
        schemaName = "gene_expression",
        queryName = "gene_expression_analysis",
        toDelete = currGEA
      )
      if (deleteGEA$rowsAffected != nrow(currGEA)) {
        stop("currGEA not deleted correctly")
      }
    }

    log_message("Writing ", nrow(de_runs), " rows to gene_expression_analysis table...")
    # push newGEA b/c listings may be different in terms of idx than old
    doneGEA <- labkey.importRows(
      baseUrl = con$config$labkey.url.base,
      folderPath = con$config$labkey.url.path,
      schemaName = "gene_expression",
      queryName = "gene_expression_analysis",
      toImport = de_runs
    )

    if (doneGEA$rowsAffected != nrow(de_runs)) {
      stop("newGEA not imported correctly")
    }

    # GEAR gets deleted and then new rows imported because
    # new mappings will be different and do not want to have leftovers
    if (nrow(de_results) > 0) {
      currGEAR <- ImmuneSpaceR:::.getLKtbl(
        con = con,
        schema = "gene_expression",
        query = "gene_expression_analysis_results"
      )

      if (nrow(currGEAR) != 0) {
        log_message("Deleting ", nrow(currGEAR), " rows from gene_expression_analysis_results...")
        delGEAR <- labkey.deleteRows(
          baseUrl = con$config$labkey.url.base,
          folderPath = con$config$labkey.url.path,
          schemaName = "gene_expression",
          queryName = "gene_expression_analysis",
          toDelete = currGEAR
        )

        postDeleteGEAR <- ImmuneSpaceR:::.getLKtbl(
          con = con,
          schema = "gene_expression",
          query = "gene_expression_analysis_results"
        )

        if (nrow(postDeleteGEAR) != 0) {
          stop("not all GEAR deleted correctly")
        }
      }

      # Import new GEAR
      de_results[is.na(de_results)] <- ""
      toImport <- data.frame(de_results, stringsAsFactors = F)

      # Load in chunks of 50000 rows
      log_message("importing ", nrow(toImport), " rows to gene_expression_analysis_results...")

      # Use labkey.query.import for importing large number of rows.
      # Depends on Rlabkey >= v2.7.0
      resGEAR <- labkey.query.import(
        baseUrl = con$config$labkey.url.base,
        folderPath = con$config$labkey.url.path,
        schemaName = "gene_expression",
        queryName = "gene_expression_analysis_results",
        toImport = toImport
      )
    }
  }
  invisible(con)
}

get_de_compatible_runs <- function(con) {
  ge_samples <- data.table(labkey.selectRows(
    baseUrl = con$config$labkey.url.base,
    folderPath = con$config$labkey.url.path,
    schemaName = "assay.expressionMatrix.matrix",
    queryName = "inputSamples_computed",
    colNameOpt = "rname",
    colSelect = c(
      "Biosample/biosample_accession",
      "Biosample/ParticipantId",
      "Biosample/arm_name",
      "Biosample/study_time_collected",
      "Biosample/study_time_collected_unit",
      "Biosample/study_accession"
    ),
    showHidden = TRUE
  ))

  # 1. Remove all arm_name * study_time_collected with less than 4 samples
  # otherwise predictive modeling cannot work
  ge_samples[,
    sample_count := unique(length(biosample_participantid)),
    by = .(biosample_arm_name, biosample_study_time_collected)
  ]
  ge_samples <- ge_samples[sample_count > 3]

  # 2. Check for baseline within each arm_name and then filter out baseline
  ge_samples[,
    includes_baseline := any(biosample_study_time_collected <= 0),
    by = .(biosample_arm_name)
  ]
  ge_samples <- ge_samples[includes_baseline == TRUE]
  ge_samples <- ge_samples[biosample_study_time_collected > 0]

  # 4. Summarize by arm_name * study_time_collected for number of subs and key
  de_runs <- unique(ge_samples[
    ,
    .(
      study_accession = biosample_study_accession,
      arm_name = biosample_arm_name,
      study_time_collected = biosample_study_time_collected,
      study_time_collected_unit = biosample_study_time_collected_unit,
      sample_count
    )
  ])
  de_runs
}
