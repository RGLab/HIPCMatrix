#' Subset Illumina raw data
#'
#' Subset Illumina expression matrix to only needed columns
#'
#' @param raw_ge_dt data.table with raw values from illumina output
#' for one sample
.subset_raw_illumina_dt <- function(raw_ge_dt){
  drop_terms <- c("bead", "array", "min", "max",
                "norm", "search", "gene", "target_id_type",
                "definition", "chromosome", "synonyms", "symbol",
                "probeid", "V([2-9]|\\d{2,3})") # allow V1
  ge_dt <- raw_ge_dt[ , grep(paste(drop_terms, collapse = "|"),
                   colnames(raw_ge_dt),
                   ignore.case = TRUE,
                   invert = TRUE),
            with = FALSE ]
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
.prep_illumina_headers <- function(raw_illumina_dt){

  detection_indices <- grep("Detection", colnames(raw_illumina_dt), ignore.case = TRUE)
  detection_names <- colnames(raw_illumina_dt)[detection_indices]
  signal_indices <- grep("Detection|ID_REF|PROBE_ID|TARGET_ID|GENE_SYMBOL",
                 colnames(raw_illumina_dt), invert = TRUE, ignore.case = TRUE)
  signal_names <- colnames(raw_illumina_dt)[signal_indices]


  change_detection <- ( all(detection_names == "Detection Pval") | all(grepl("P_VALUE|PVAL", detection_names)) ) &
    length(detection_names) > 0
  change_signal <- !all( grepl("AVG_Signal", signal_names) ) & !all( grepl("^ES\\d{6,7}$", signal_names) )
  change_probeid <- !any( grepl("ID_REF", colnames(raw_illumina_dt)) )

  if ( change_detection ) {

    if ( all(detection_names == "Detection Pval") ) {
      detection_names <- paste0(signal_names, ".", detection_names)
    }
    detection_names <- gsub("DETECTION_(PVAL|P_VALUE)", "Detection Pval", detection_names)
    setnames(raw_illumina_dt, detection_indices, detection_names)

  }

  if ( change_signal ) {
    if ( any(grepl("RAW", signal_names)) ) {
      signal_names <- gsub("RAW_SIGNAL", "AVG_Signal", signal_names)
    } else if (!any(grepl("SAMPLE", signal_names))) {
      signal_names <- paste0(signal_names, ".AVG_Signal")
    } else if (all(grepl("^SAMPLE", signal_names))){
      signal_names <- paste0(signal_names, ".AVG_Signal")
    } else {
      signal_names <- gsub("AVG_Signal", "SAMPLE", signal_names)
    }
    signal_names <- gsub("SIGNAL", "Signal", signal_names)
    setnames(raw_illumina_dt, signal_indices, signal_names)
  }

  if ( change_probeid ) {
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
#' @param rawFl file path to raw illumina file
#'
.processIllumina <- function(rawFl){
  em <- fread(rawFl)

  # check for known issues that would hinder background correction
  # 1. Subjects with no fluorescence measurements
  badSubs <- apply(em, 2, function(x){ all( x == 0 ) })
  if (any(badSubs)) {
    nms <- names(badSubs[ badSubs == TRUE ])
    es <- regmatches(nms, regexpr("(ES|GSM)\\d{6,7}", nms))
    em <- em[ , grep(es, colnames(em), invert = TRUE), with = FALSE]
  }

  # 2. Control or misnamed probes (not unique) - e.g. "NEGATIVE"
  em <- em[ grep("ILMN", em$ID_REF) ]
  write.table(em, rawFl, sep = "\t", row.names = FALSE, quote = FALSE)

  # Can only background correct using detection pvals.
  # Immport-derived files may already have this done in some cases.
  if (any(grepl("Detection", colnames(em)))){
    esList <- limma::read.ilmn(rawFl,
                               expr = "AVG_Signal",
                               probeid = "ID_REF")
    em <- data.table(limma::nec(esList)$E, keep.rownames = TRUE)
  }

  # Fix names for future mapping if necessary as read.ilmn()
  # leaves a suffix on `<smpl>_AVG_Signal` to be `<smpl>_`
  tags <- "BS|GSM|ES"
  if (any(grepl(tags, colnames(em)))) {
    nmsVals <- grep(tags, colnames(em), value = TRUE)
    rep  <- gsub("_", "", nmsVals) # SDY162
    setnames(em, nmsVals, rep)
  }

  return(em)
}
