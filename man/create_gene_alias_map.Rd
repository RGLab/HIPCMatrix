% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_symbol_mapping.R
\name{create_gene_alias_map}
\alias{create_gene_alias_map}
\title{Create Gene Alias Map}
\usage{
create_gene_alias_map()
}
\value{
data.table with four columns: \describe{
\item{SYMBOL}{"official" gene symbol}
\item{ALIAS}{Aliases for the official gene symbol. Each entry in this
column should be unique}
\item{HGNC}{HGNC ID}
\item{ENTREZ}{Entrez ID}
}
}
\description{
This will download the current gene-alias mapping info from
the HUGO gene nomenclature committee and create a mapping from
all recorded gene aliases to official gene symbols.
}
\details{
Some aliases map to multiple gene symbols. There are a couple of
different cases, which are dealt with differently:
\enumerate{
\item A gene alias maps to multiple gene symbols including itself:
Remove mappings to other symbols.
\item A gene alias maps to multiple gene symbols not including itself:
Remove this alias.
}
}
\references{
\href{https://www.genenames.org/}{HGNC}
}
