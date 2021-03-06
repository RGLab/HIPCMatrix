% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_files.R
\name{.map_sampleid_to_biosample_accession}
\alias{.map_sampleid_to_biosample_accession}
\title{Map experiment-sample or geo accessions to biosample accessions}
\usage{
.map_sampleid_to_biosample_accession(exprs_dt, gef)
}
\arguments{
\item{exprs_dt}{data.table of gene expression with one column per sample,
one row per feature.}

\item{gef}{result of ISCon$getDataset("gene_expression_files") for one run.}
}
\description{
Map experiment-sample or geo accessions to biosample accessions
}
