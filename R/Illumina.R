#' Subset Illumina EM
#'
#' Subset Illumina expression matrix to only needed columns
#'
#' @param em expression
.subsetIlluminaEM <- function(em){
  badTerms <- c("bead", "array", "min", "max",
                "norm", "search", "gene", "target_id_type",
                "definition", "chromosome", "synonyms", "symbol",
                "probeid", "V([2-9]|\\d{2,3})") # allow V1
  em <- em[ , grep(paste(badTerms, collapse = "|"), colnames(em), ignore.case = TRUE, invert = TRUE), with = FALSE ]
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
#' @param em expression
.prepIlluminaHeaders <- function(em){
  detLoc <- grep("Detection", colnames(em), ignore.case = TRUE)
  detVals <- colnames(em)[detLoc]
  nmsLoc <- grep("Detection|ID_REF|PROBE_ID|TARGET_ID|GENE_SYMBOL",
                 colnames(em), invert = TRUE, ignore.case = TRUE)
  nmsVals <- colnames(em)[nmsLoc]
  chgDet <- (all(detVals == "Detection Pval") | all(grepl("P_VALUE|PVAL", detVals))) &
    length(detVals) > 0
  chgNms <- !all(grepl("AVG_Signal", nmsVals)) & !all(grepl("^ES\\d{6,7}$", nmsVals))
  chgPrb <- !any(grepl("ID_REF", colnames(em)))

  if (chgDet) {
    if (all(detVals == "Detection Pval")) {
      detVals <- paste0(nmsVals, ".", detVals)
    }
    detVals <- gsub("DETECTION_(PVAL|P_VALUE)", "Detection Pval", detVals)
    setnames(em, detLoc, detVals)
  }

  if (chgNms) {
    if ( any(grepl("RAW", nmsVals))) {
      nmsVals <- gsub("RAW_SIGNAL", "AVG_Signal", nmsVals)
    } else if (!any(grepl("SAMPLE", nmsVals))) {
      nmsVals <- paste0(nmsVals, ".AVG_Signal")
    } else if (all(grepl("^SAMPLE", nmsVals))){
      nmsVals <- paste0(nmsVals, ".AVG_Signal")
    } else {
      nmsVals <- gsub("AVG_Signal", "SAMPLE", nmsVals)
    }
    nmsVals <- gsub("SIGNAL", "Signal", nmsVals)
    setnames(em, nmsLoc, nmsVals)
  }

  if (chgPrb) {
    prb <- grep("PROBE_ID|V1|TARGET_ID", colnames(em))
    setnames(em, prb, "ID_REF")
  }

  return(em)
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
