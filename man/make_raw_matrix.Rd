% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_raw_matrix.R
\name{make_raw_matrix}
\alias{make_raw_matrix}
\title{Make Raw Matrix}
\usage{
make_raw_matrix(platform, gef, input_files, verbose = FALSE)
}
\arguments{
\item{platform}{Illumina, Affymetrix, or NA}

\item{gef}{result of \code{ISCon$getDataset("gene_expression_files")} for one run.}

\item{input_files}{input file names}

\item{verbose}{print verbose logging statements?}
}
\value{
A data.table containing raw counts or background-corrected probe
intensities with a feature_id column and one column per biosample_accession
}
\description{
Using raw input files, create a matrix with raw counts or background-corrected
probe intensities, mapped to biosample accessions.
}
