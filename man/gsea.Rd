% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_set_enrichment_analysis.R
\name{gsea}
\alias{gsea}
\title{Gene Set Enrichment Analysis}
\usage{
gsea(
  eset,
  set_name = "msigdb",
  gene_sets = NULL,
  contrast = "study_time_collected",
  baseline = NULL
)
}
\arguments{
\item{eset}{expressionset}

\item{set_name}{Name of predefined set of gene signatures. Choose from:
\code{chaussabel}, \code{blood_transcription}, \code{msigdb}}

\item{gene_sets}{A list of vectors of gene names, each entry corresponding
to a gene set. If specified, this will be used in place of the "set_name"
argument to test gene sets.}

\item{contrast}{contrast}

\item{baseline}{Value to compare other contrast values against. By default,
uses first when coerced to factor.}
}
\description{
Perform a gene set enrichment analysis on an expressionset using the
CAMERA method.
}
