#' data.table mapping gene aliases to gene symbols
#'
#' Gene-alias mapping info from
#' the HUGO gene nomenclature committee.
#'
#' @format data.table with four columns: \describe{
#' \item{SYMBOL}{"official" gene symbol}
#' \item{ALIAS}{Aliases for the official gene symbol. Each entry in this
#' column should be unique}
#' \item{HGNC}{HGNC ID}
#' \item{ENTREZ}{Entrez ID}
#' }
#'
#' @details Some aliases map to multiple gene symbols. There are a couple of
#' different cases, which are dealt with differently:
#' \enumerate{
#'   \item A gene alias maps to multiple gene symbols including itself:
#'   Remove mappings to other symbols.
#'   \item A gene alias maps to multiple gene symbols not including itself:
#'   Remove this alias.
#' }
#'
#' @references \href{https://www.genenames.org/}{HGNC}
"hgncAlias2Symbol"

#' Date hgncAlias2Symbol update
"hgncAlias2Symbol_version"


download_hgnc_complete_set <- function(file = tempfile()) {
  url <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json"
  utils::download.file(url, file)
  tmp <- jsonlite::fromJSON(file)
  data <- data.table(tmp$response$docs)
  return(data)
}

create_gene_alias_map_from_hgnc_set <- function(data) {
  data <- data[, list(symbol, alias_symbol, prev_symbol, hgnc_id, entrez_id)]
  # add rows for all previous and current aliases
  data <- apply(data, 1, function(x) {
    aliases <- c(unlist(x$alias_symbol), unlist(x$prev_symbol), x$symbol)
    newRows <- data.table(
      SYMBOL = x$symbol,
      ALIAS = aliases,
      HGNC = as.numeric(gsub("HGNC:", "", x$hgnc_id)),
      ENTREZ = x$entrez_id
    )
  })
  data <- rbindlist(data)
  setorder(data, "SYMBOL")
  data <- data[!duplicated(data)]


  # Handle genes with entries both as an alias and as a symbol
  # by removing mappings to other mappings and maintaining the 'symbol' observation
  # based on the premise that mappings to other symbols are a historic artifact.
  # This represents approximately 0.7% of all aliases as of 12/2020
  selfMap <- data[ALIAS == SYMBOL]
  otherMap <- data[ALIAS != SYMBOL]
  overlap <- intersect(selfMap$ALIAS, otherMap$ALIAS)
  toRemove <- which((data$ALIAS %in% overlap) & (data$ALIAS != data$SYMBOL))
  data <- data[!toRemove]

  # Some aliases still map to multiple symbols
  # Remove these aliases since we have no good way to
  # tell which symbol is the most accurate.
  # This represents approximately 1.3% of all aliases as of 12/2020
  multiMapped <- data[, .N, ALIAS][N > 1]
  data <- data[!ALIAS %in% multiMapped$ALIAS]

  if (length(unique(data$ALIAS)) != nrow(data)) {
    stop("There are multi-mappings of alias to symbol. Must correct and re-run!")
  }

  data
}

#' Create Gene Alias Map
#'
#' This will download the current gene-alias mapping info from
#' the HUGO gene nomenclature committee and create a mapping from
#' all recorded gene aliases to official gene symbols.
#'
#' @return data.table with four columns: \describe{
#' \item{SYMBOL}{"official" gene symbol}
#' \item{ALIAS}{Aliases for the official gene symbol. Each entry in this
#' column should be unique}
#' \item{HGNC}{HGNC ID}
#' \item{ENTREZ}{Entrez ID}
#' }
#'
#' @details Some aliases map to multiple gene symbols. There are a couple of
#' different cases, which are dealt with differently:
#' \enumerate{
#'   \item A gene alias maps to multiple gene symbols including itself:
#'   Remove mappings to other symbols.
#'   \item A gene alias maps to multiple gene symbols not including itself:
#'   Remove this alias.
#' }
#'
#' @references \href{https://www.genenames.org/}{HGNC}
#' @export
create_gene_alias_map <- function() {
  data <- download_hgnc_complete_set()
  create_gene_alias_map_from_hgnc_set(data)
}
