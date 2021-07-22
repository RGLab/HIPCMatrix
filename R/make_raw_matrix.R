#' Make Raw Matrix
#'
#' Using raw input files, create a matrix with raw counts or background-corrected
#' probe intensities, mapped to biosample accessions.
#'
#' @return A data.table containing raw counts or background-corrected probe
#' intensities with a feature_id column and one column per biosample_accession
#'
#' @param platform Illumina, Affymetrix, or NA
#' @param gef result of \code{ISCon$getDataset("gene_expression_files")} for one run.
#' @param input_files input file names
#' @param verbose print verbose logging statements?
#' @export
make_raw_matrix <- function(platform,
                            gef,
                            input_files,
                            verbose = FALSE) {


  if (verbose) log_message("Creating matrix of raw expression...")
  # Generate background corrected raw matrices for affy and illumina
  # For RNAseq pass through raw counts file.
  exprs_dt <- switch(
    platform,
    "Affymetrix" = .process_affy(
      input_files,
      verbose = verbose
    ),
    "Illumina" = .process_illumina(
      input_files,
      verbose = verbose
    ),
    "NA" = .process_rna_seq(
      input_files,
      verbose = verbose
    ),
    fread(input_files)
  )



  # Ensure probe column is named 'feature_id'
  exprs_dt <- .map_feature_id_col(exprs_dt)

  # Ensure all probes have names. Not a problem for most studies.
  # Note that values can still be NA here and may be due to a handful
  # of samples having problems (e.g. SDY224 / SDY212).  These NA values
  # are removed during normalization.
  exprs_dt <- exprs_dt[!is.na(feature_id) & feature_id != ""]

  # Map expsample or geo accession to biosample accession
  exprs_dt <- .map_sampleid_to_biosample_accession(exprs_dt, gef)

  # Subset to biosamples in selected_biosamples from UI
  exprs_dt <- exprs_dt[, colnames(exprs_dt) %in% c("feature_id", gef$biosample_accession),
    with = FALSE
  ]

  # Check that all gef samples are in the matrix - may not be the case if some
  # were removed for QC issues.  User should re-run matrix with these samples
  # removed.
  if (!all(gef$biosample_accession %in% colnames(exprs_dt))) {
    stop(
      "Some selected biosamples are not found in raw-matrix.",
      "These may have been removed for QC issues. Please check and re-run"
    )
  }

  # Ensure colOrder
  sampleids <- grep("feature_id", names(exprs_dt), invert = TRUE)
  feature_id <- grep("feature_id", names(exprs_dt))
  setcolorder(exprs_dt, c(feature_id, sampleids))

  return(exprs_dt)
}


# ----- Illumina Helpers -----


#' Subset Illumina raw data
#'
#' Subset Illumina expression matrix to only needed columns
#'
#' @param raw_ge_dt data.table with raw values from illumina output
#' for one sample
.subset_raw_illumina_dt <- function(raw_ge_dt) {
  drop_terms <- c(
    "bead", "array", "min", "max",
    "norm", "search", "gene", "target_id_type",
    "definition", "chromosome", "synonyms", "symbol",
    "probeid", "V([2-9]|\\d{2,3})"
  ) # allow V1
  ge_dt <- raw_ge_dt[,
    grep(paste(drop_terms, collapse = "|"),
      colnames(raw_ge_dt),
      ignore.case = TRUE,
      invert = TRUE
    ),
    with = FALSE
  ]
  ge_dt
}

#' Prep Illumina Headers
#'
#' For read.ilmn() to work correctly the signal cols need format <smpl>.<exprValTerm>
#' with the expression-value term being something like "AVG_Signal" or "SAMPLE". Then the
#' detection p-value cols must have format <smpl>.Detection Pval. This way the
#' rawElist creates the $E and $other$Detection matrices with same colnames of <smpl>.
#' Since in this script the read.ilmn(file, exprs = "AVG_Signal", probeid = "ID_REF) is
#' hardcoded. The vars are substituted here for other viable versions, e.g. SAMPLE and PROBE_ID.
#'
#' @param raw_illumina_dt data.table of raw illumina values.
#'
#' @return data.table with cleaned up names ready for \code{read.ilmn()}
.prep_illumina_headers <- function(raw_illumina_dt) {
  detection_indices <- grep("Detection", colnames(raw_illumina_dt), ignore.case = TRUE)
  detection_names <- colnames(raw_illumina_dt)[detection_indices]
  signal_indices <- grep("Detection|ID_REF|PROBE_ID|TARGET_ID|GENE_SYMBOL",
    colnames(raw_illumina_dt),
    invert = TRUE, ignore.case = TRUE
  )
  signal_names <- colnames(raw_illumina_dt)[signal_indices]


  change_detection <- (all(detection_names == "Detection Pval") | all(grepl("P_VALUE|PVAL", detection_names))) &
    length(detection_names) > 0
  change_signal <- !all(grepl("AVG_Signal", signal_names)) & !all(grepl("^ES\\d{6,7}$", signal_names))
  change_probeid <- !any(grepl("ID_REF", colnames(raw_illumina_dt)))

  if (change_detection) {
    if (all(detection_names == "Detection Pval")) {
      detection_names <- paste0(signal_names, ".", detection_names)
    }
    detection_names <- gsub("DETECTION_(PVAL|P_VALUE)", "Detection Pval", detection_names)
    setnames(raw_illumina_dt, detection_indices, detection_names)
  }

  if (change_signal) {
    if (any(grepl("RAW", signal_names))) {
      signal_names <- gsub("RAW_SIGNAL", "AVG_Signal", signal_names)
    } else if (!any(grepl("SAMPLE", signal_names))) {
      signal_names <- paste0(signal_names, ".AVG_Signal")
    } else if (all(grepl("^SAMPLE", signal_names))) {
      signal_names <- paste0(signal_names, ".AVG_Signal")
    } else {
      signal_names <- gsub("AVG_Signal", "SAMPLE", signal_names)
    }
    signal_names <- gsub("SIGNAL", "Signal", signal_names)
    setnames(raw_illumina_dt, signal_indices, signal_names)
  }

  if (change_probeid) {
    probeid_index <- grep("PROBE_ID|V1|TARGET_ID", colnames(raw_illumina_dt))
    setnames(raw_illumina_dt, probeid_index, "ID_REF")
  }

  return(raw_illumina_dt)
}

#' Process Illumina
#'
#' For Illumina, correct background using nec() which
#' calls normexp.fit.control and normexp.signal to use
#' negative controls as identified by detection p-values
#' to remove noise.
#'
#' @param raw_file_path file path to raw illumina file
#' @param verbose print verbose logging statements?
#'
.process_illumina <- function(raw_file_path,
                              verbose = FALSE) {
  if (verbose) log_message("Processing illumina files...")
  raw_dt <- fread(raw_file_path)

  # check for known issues that would hinder background correction
  # 1. Subjects with no fluorescence measurements
  badSubs <- apply(raw_dt, 2, function(x) {
    all(x == 0)
  })
  if (any(badSubs)) {
    nms <- names(badSubs[badSubs == TRUE])
    es <- regmatches(nms, regexpr("(ES|GSM)\\d{6,7}", nms))
    raw_dt <- raw_dt[, grep(es, colnames(raw_dt), invert = TRUE), with = FALSE]
  }

  # 2. Control or misnamed probes (not unique) - e.g. "NEGATIVE"
  raw_dt <- raw_dt[grep("ILMN", raw_dt$ID_REF)]
  utils::write.table(raw_dt, raw_file_path, sep = "\t", row.names = FALSE, quote = FALSE)

  # Can only background correct using detection pvals.
  # Immport-derived files may already have this done in some cases.
  if (any(grepl("Detection", colnames(raw_dt)))) {
    # Get intensities
    if (verbose) log_message("Reading in background-corrected probe intensities...")
    esList <- limma::read.ilmn(raw_file_path,
      expr = "AVG_Signal",
      probeid = "ID_REF",
      verbose = FALSE
    )
    # Background correction
    raw_dt <- data.table(limma::nec(esList)$E, keep.rownames = TRUE)
  }

  # Fix names for future mapping if necessary as read.ilmn()
  # leaves a suffix on `<smpl>_AVG_Signal` to be `<smpl>_`
  tags <- "BS|GSM|ES"
  if (any(grepl(tags, colnames(raw_dt)))) {
    nmsVals <- grep(tags, colnames(raw_dt), value = TRUE)
    rep <- gsub("_", "", nmsVals) # SDY162
    setnames(raw_dt, nmsVals, rep)
  }

  return(raw_dt)
}


# ----- Affymetrix -----

#' Process Affy
#'
#' @param input_files input file names
#' @param verbose print verbose logging statements?
#'
#' @return background-corrected data.table of probe intensities
#'
#' @details
#' Requires annotation packages which can be installed from bioconductor.
#' \code{affy::justRMA()} will use the appropriate package.
#' Two are customCDF packages loaded from UpdateAnno: \code{huex10stv2cdf}
#' and \code{hursta2a520709cdf}
#'
.process_affy <- function(input_files,
                          verbose = FALSE) {
  if (verbose) log_message("Processing affymetrix files...")
  if (verbose) log_message("Processing ", length(input_files), " CEL files")
  # Background Correction Notes:
  # 'background' = TRUE performs function similar to normexp.fit.control and normexp.signal
  # from limma package.
  wd <- getwd()
  setwd("/") # b/c filepaths are absolute and justRMA prepends wd
  eset <- try(affy::justRMA(filenames = input_files, normalize = FALSE, background = TRUE),
    silent = TRUE
  )
  setwd(wd)
  if ("try-error" %in% class(eset)) stop(eset)

  exprs_dt <- data.table(Biobase::exprs(eset), keep.rownames = TRUE)
  setnames(exprs_dt, "rn", "feature_id")

  # Names come from input_files. In case of isGeo, these are not exact
  # matches but usually have the gsm accession in them.
  if (any(grep("GSM", colnames(exprs_dt)))) {
    sample_names <- grep(
      "feature_id",
      colnames(exprs_dt),
      invert = TRUE,
      value = TRUE
    )
    gsms <- regmatches(sample_names, regexpr("GSM\\d{6,7}", sample_names))
    setnames(exprs_dt, sample_names, gsms)
  }

  return(exprs_dt)
}

#' Process RNA-seq
#'
#' Somewhat redundant in case of GEO files as there is no actual processing
#' of raw counts files ...
#'
#' @param input_files input file names
#' @param verbose print verbose logging statements?
.process_rna_seq <- function(input_files,
                             verbose = FALSE) {
  if (verbose) log_message("Processing RNA-seq files...")
  lf <- lapply(input_files, fread)
  exprs <- data.table(Reduce(f = function(x, y) {
    merge(x, y, all = TRUE)
  }, lf))
}


#' Standardize probe column name
#'
#' Different platforms use different names for the feature id column.
#' This will update this column to be called \code{feature_id}
#'
#' @param exprs_dt data.table of gene expression with one column per sample,
#' one row per feature.
#'
.map_feature_id_col <- function(exprs_dt) {
  if (!any(grepl("feature_id", colnames(exprs_dt)))) {

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

#' Process Two Color Array
#'
#' Two color array processing using limma and assuming genepix files.
#' bc.method is the background correction method and normexp is used to match
#' work with Illumina and Affymetrix.
#'
#' @param path path to input files
#' @param verbose print verbose logging statements?
#'
#' @export
.process_two_color_array <- function(path,
                                     verbose = FALSE) {
  if (verbose) log_message("Processing two-color-array...")
  RG <- limma::read.maimages(files = path, source = "genepix")
  MA <- limma::normalizeWithinArrays(RG, bc.method = "normexp", method = "none")
  ge_dt <- data.table(feature_id = MA$genes$ID, gsm = MA$A[, 1])

  # RM dup probes due to multiple probes per spot
  # RM NA vals possibly from background correction issues
  ge_dt <- ge_dt[!duplicated(feature_id) & !is.na(gsm)]

  # RM unmappable feature_ids
  ge_dt <- ge_dt[grep("EMPTY|bsid", feature_id, invert = TRUE)]
}
