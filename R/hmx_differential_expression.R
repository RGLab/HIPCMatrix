#' @include HMX.R
NULL


#' Compute Differential Expression using Emperical Bayes Statistics
#'
#' Compute differential expression over time on microarray data using the
#' \code{\link[limma]{eBayes}} method from the \code{limma} package. It
#' reports all differentially expressed genes for each timepoint.
#'
#' @return list with one entry for each
#'
#' @param eset expressionset object with normalized microarray expression.
#' @param p_val_cutoff adjusted p-value cutoff for determining differentially
#' expressed genes.
#'
#' @export
find_de_genes_eBayes <- function(
                                 eset,
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
    message(paste0(run, " has fewer than 2 participants with multiple timepoints! Skipping."))
    return()
  }

  mm <- model.matrix(formula("~participant_id + coef"), pd)

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
    message(paste0("Linear model not able to be fit for ", run, ". Skipping to next matrix"))
    return()
  }

  coefs <- grep("^coef", colnames(mm), value = TRUE)
  tt_list <- list()
  for (coef in coefs) {

    # Check that coef can be used
    tt <- data.table(limma::topTable(fit, coef = coef, number = Inf))
    if (all(is.na(tt$adj.P.Val))) {
      message(paste0(coef, " has all NA values for adj.P.Val. Skipping to next coef."))
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

#' runGEAnalysis
#'
#' @description Run Differential Expression analysis on matrices in a study
#'
#' @param rerun Force re-run if results already in cache?
#'
HMX$set(
  which = "public",
  name = "runGEAnalysis",
  value = function(rerun = FALSE) {
    if (grepl("^IS", self$study)) {
      stop("runGEAnalysis is not designed for ImmuneSignatures!")
    }

    if (!grepl("SDY", self$study)) {
      stop("run runGEAnalysis one study at a time!")
    }

    if ("de_results" %in% names(self$cache) & !rerun) {
      return(self$cache$de_results)
    }

    if (self$config$verbose) {
      log_message("Running Differential Expression Analysis...")
    }

    self$getGEInputs()
    contrast <- c("study_time_collected", "study_time_collected_unit")
    coefs <- unique(self$cache$GE_inputs[, c("arm_name", contrast), with = FALSE])

    if (!any(coefs$study_time_collected <= 0)) {
      log_message("No baseline timepoints available in any cohort. Analysis not run.")
      return()
    }

    if (sum(coefs$study_time_collected > 0) == 0) {
      log_message("No post-baseline timepoints available in any cohort. Analysis not run.")
      return()
    }

    GEA_list <- vector("list")
    GEAR_list <- vector("list")

    runs <- self$cache$GE_matrices$name
    gef <- self$getDataset("gene_expression_files", original_view = T)

    idx <- 1 # analysis accession key
    for (run in runs) {
      eset <- self$getGEMatrix(matrixName = run, outputType = "normalized", annotation = "latest")
      tt_list <- find_de_genes_eBayes(eset)


      cm <- self$getDataset("cohort_membership")
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

    self$cache[["de_runs"]] <- de_runs
    self$cache[["de_results"]] <- de_results
    return(de_results)
  }
)

#' uploadGEAnalysisResults
#'
#' @description Upload differential gene expression analysis results to
#' server. This updates the gene_expression_analysis and the
#' gene_expression_analysis_results tables with updated results.
HMX$set(
  which = "public",
  name = "uploadGEAnalysisResults",
  value = function() {
    if (!grepl("SDY", self$study)) {
      stop("Please run for one study at a time.")
    }
    de_results <- self$runGEAnalysis()
    de_runs <- self$cache$de_runs
    if (nrow(de_results) == 0) {
      stop("No differentially expressed genes found.")
    }

    if (nrow(de_runs) > 0) {

      # delete old GEA
      currGEA <- tryCatch(
        ImmuneSpaceR:::.getLKtbl(
          con = self,
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
          baseUrl = self$config$labkey.url.base,
          folderPath = self$config$labkey.url.path,
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
        baseUrl = self$config$labkey.url.base,
        folderPath = self$config$labkey.url.path,
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
          con = self,
          schema = "gene_expression",
          query = "gene_expression_analysis_results"
        )

        if (nrow(currGEAR) != 0) {
          log_message("Deleting ", nrow(currGEAR), " rows from gene_expression_analysis_results...")
          delGEAR <- labkey.deleteRows(
            baseUrl = self$config$labkey.url.base,
            folderPath = self$config$labkey.url.path,
            schemaName = "gene_expression",
            queryName = "gene_expression_analysis",
            toDelete = currGEAR
          )

          postDeleteGEAR <- ImmuneSpaceR:::.getLKtbl(
            con = self,
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
        log_message("importing ", nrow(toImport), " rows to gene_expression_analysis_restuls...")

        # Use labkey.query.import for importing large number of rows.
        # Depends on Rlabkey >= v2.7.0
        resGEAR <- labkey.query.import(
          baseUrl = self$config$labkey.url.base,
          folderPath = self$config$labkey.url.path,
          schemaName = "gene_expression",
          queryName = "gene_expression_analysis_results",
          toImport = toImport
        )
      }
    }
    invisible(self)
  }
)

HMX$set(
  which = "public",
  name = "checkImpliedGEAR",
  value = function() {
    impliedGEA <- data.table(labkey.selectRows(
      baseUrl = self$config$labkey.url.base,
      folderPath = self$config$labkey.url.path,
      schemaName = "assay.expressionMatrix.matrix",
      queryName = "inputSamples",
      colNameOpt = "rname",
      showHidden = TRUE
    ))

    # 1. Remove all arm_name * study_time_collected with less than 4 replicates
    # otherwise predictive modeling cannot work
    impliedGEA[, subs := unique(length(biosample_participantid)), by = .(biosample_arm_name, biosample_study_time_collected)]
    impliedGEA <- impliedGEA[subs > 3]

    # 2. Check for baseline within each arm_name and then filter out baseline
    impliedGEA[, baseline := any(biosample_study_time_collected <= 0), by = .(biosample_arm_name)]
    impliedGEA <- impliedGEA[baseline == TRUE]
    impliedGEA <- impliedGEA[biosample_study_time_collected > 0]

    # 3. Generate key
    impliedGEA[, key := paste(biosample_arm_name, biosample_study_time_collected, biosample_study_time_collected_unit)]

    # 4. Summarize by arm_name * study_time_collected for number of subs and key
    smryGEA <- impliedGEA[, list(key = unique(key), subs = unique(subs)), by = .(biosample_arm_name, biosample_study_time_collected)]
    dim(smryGEA)[[1]] > 0
  }
)
