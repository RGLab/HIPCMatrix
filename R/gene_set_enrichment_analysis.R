#' Gene Set Enrichment Analysis
#'
#' Perform a gene set enrichment analysis on an expressionset using the
#' CAMERA method.
#'
#' @param eset expressionset
#' @param set_name Name of predefined set of gene signatures. Choose from:
#' \code{chaussabel}, \code{blood_transcription}, \code{msigdb}
#' @param gene_sets  A list of vectors of gene names, each entry corresponding
#' to a gene set. If specified, this will be used in place of the "set_name"
#' argument to test gene sets.
#' @param contrast contrast
#' @param baseline Value to compare other contrast values against. By default,
#' uses first when coerced to factor.
#' @export
gsea <- function(eset,
                 set_name = "msigdb",
                 gene_sets = NULL,
                 contrast = "study_time_collected",
                 baseline = NULL) {

  if (is.null(gene_sets)) {
    gene_sets <- switch(
      set_name,
      "chaussabel" = chaussabel_modules,
      "blood_transcription" = emory_blood_transcript_modules,
      "msigdb" = msigdb_immunologic_signatures,
      stop("Invalid signature set.")
    )
  }

  indices <- limma::ids2indices(gene_sets, Biobase::fData(eset)$gene_symbol)
  Biobase::pData(eset)[, contrast] <- as.factor(Biobase::pData(eset)[, contrast])

  if (length(unique(Biobase::pData(eset)[, contrast])) < 2) {
    stop("Cohort does not have multiple timepoints for differential expression data. Unable to analyze.")
  }

  # Relevel contrast if baseline is not null
  if (!is.null(baseline)) {
    pd <- Biobase::pData(eset)
    if (!baseline %in% pd[, contrast]) stop(baseline, " not found in ", contrast)
    values <- pd[, contrast]
    values <- factor(values, levels = c(baseline, unique(as.character(values[values != baseline]))))
    pd[, contrast] <- values
    Biobase::pData(eset) <- pd
  }

  mm <- stats::model.matrix(
    stats::formula(paste0("~participant_id + ", contrast)),
    eset
  )
  colnames(mm) <- make.names(colnames(mm))



  lev <- levels(Biobase::pData(eset)[, contrast])
  lev_label <- paste0(contrast, lev[2:length(lev)])
  lev <- make.names(lev_label)
  cam_list <- vector("list", length = length(lev))

  for (i in 1:length(lev)) {
    contrasts <- limma::makeContrasts(contrasts = lev[[i]], levels = mm)
    res <- limma::camera(
      eset,
      index = indices,
      design = mm,
      contrast = contrasts,
      allow.neg.cor = TRUE,
      inter.gene.cor = NULL
    )

    dt <- data.table(res)
    dt[, Module := rownames(res)]
    dt[, Coefficient := lev_label[[i]]]
    dt[, PValue_10log10 := ifelse(Direction == "Up", -10 * log10(PValue), 10 * log10(PValue))]
    cam_list[[i]] <- dt
  }


  dt <- rbindlist(cam_list)
  setcolorder(dt, c("Module", "Coefficient", "NGenes", "Correlation", "Direction", "PValue_10log10", "FDR"))
  dt <- setorder(dt, FDR)
  dt[, Coefficient := gsub(contrast, "", Coefficient)]

  dt
}


# Look for documentation in HMX$run_gsea
run_gsea <- function(con,
                     matrix_name = NULL,
                     cohort_type = NULL,
                     set_name = "msigdb",
                     gene_sets = NULL) {
  if (is.null(matrix_name) & is.null(cohort_type)) {
    stop("matrix_name or cohort_type must be specified")
  }

  if (!is.null(matrix_name) & length(matrix_name) > 1) {
    stop("matrix_name must be length 1")
  } else if (!is.null(cohort_type) & length(cohort_type) > 1) {
    stop("cohort_type must be length 1")
  }

  eset <- con$getGEMatrix(
    matrixName = matrix_name,
    cohortType = cohort_type
  )
  eset$timepoint <- paste0(eset$study_time_collected, " ", eset$study_time_collected_unit)
  eset$timepoint <- factor(eset$timepoint)

  gsea_result <- gsea(
    eset = eset,
    set_name = set_name,
    gene_sets = gene_sets,
    contrast = "timepoint"
  )
}
