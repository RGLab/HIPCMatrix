

test_that("download_hgnc_complete_set returns correct format", {
  data <- HIPCMatrix:::download_hgnc_complete_set()
  expect_s3_class(data, "data.table")
  expect_true(all(
    c("symbol",
      "alias_symbol",
      "prev_symbol",
      "hgnc_id",
      "entrez_id") %in% names(data)))
})

# Test data for creating mapping table
data <- data.table(
  symbol = letters[1:20],
  alias_symbol = rep(list(c(NULL)), 20),
  prev_symbol = rep(list(c(NULL)), 20),
  hgnc_id = paste0("HGNC:", 1:20),
  entrez_id = as.character(51:70)
)
# Add some aliases

data[symbol == "b", alias_symbol := c("b1", "b2")]
data[symbol == "c", prev_symbol := c("1c", "5c")]

# edge case: An alias maps to multiple symbols, including itself
data[symbol == "f", alias_symbol := c("f")]
data[symbol == "g", alias_symbol := c("f")]

# edge case: An alias maps to multiple symbols not includig itself
data[symbol == "j", alias_symbol := "gene"]
data[symbol == "k", alias_symbol := "gene"]

test_that("create_gene_alias_map_from_hgnc_set returns correct format", {
  map <- HIPCMatrix:::create_gene_alias_map_from_hgnc_set(data)
  expect_s3_class(map, "data.table")
  expect_equal(names(map), c("SYMBOL", "ALIAS", "HGNC", "ENTREZ"))
  expect_true(all(apply(map, 2, class) == "character"))
  expect_lt(length(unique(map$SYMBOL)), nrow(map))
  expect_equal(length(unique(map$ALIAS)), nrow(map))
})

test_that("create_gene_alias_map_from_hgnc_set handles edge cases", {
  map <- HIPCMatrix:::create_gene_alias_map_from_hgnc_set(data)

  # edge case: An alias maps to multiple symbols, including itself. use self-mapping
  expect_equal(nrow(map[ALIAS == "f"]), 1)
  expect_equal(map[ALIAS == "f", SYMBOL], "f")

  # edge case: An alias maps to multiple symbols, not including itself. Drop.
  expect_equal(nrow(map[ALIAS == "gene"]), 0)
})


