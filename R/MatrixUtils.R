
#' write Single Matrix
#'
#' @param ge_tbl object containing gene expression values (matrix, data.frame, etc)
#' @param supp_files_dir path to base directory (via \code{.get_supp_files_dir()})
#' @param study study accession eg \code{SDY269}
.write_single_ge_mx <- function(ge_tbl, supp_files_dir, study){
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
  input_files <- .write_single_ge_mx(ge_df, supp_files_dir, study)
}
#######################################
###            MAPPING              ###
#######################################

#' Standardize probe column name
#'
#' @param exprs data.table of expression
#'
.map_feature_id_col <- function(exprs){
  if (!any(grepl("feature_id", colnames(exprs)))) {

    # If Illumina from Immport
    prbCol <- grep("id_ref", colnames(exprs), ignore.case = TRUE)

    # If RNAseq then accept gene* or V1 col
    if (length(prbCol) == 0) {
      prbCol <- grep("gene|^V1$", colnames(exprs), ignore.case = TRUE)
    }

    # In case of features in rownames, e.g. from GEO
    if (length(prbCol) == 0) {
      prbCol <- "rn"
    }

    setnames(exprs, prbCol, "feature_id")
  }
  return(exprs)
}
