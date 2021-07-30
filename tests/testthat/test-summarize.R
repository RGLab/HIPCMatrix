
test_that("summarize_by_gene_symbol mean", {
  # Genes with some duplicates to mimic features that map to multiple genes
  genes <- c(NA, "gene1", "gene2", "gene1", "gene3", "gene4", "gene1", NA)
  features <- letters[1:8]
  feature_gene_map <- data.table(
    featureid = features,
    genesymbol = genes
  )

  # fake gene expression data with features in a random order
  ge_norm <- data.table(
    feature_id = sample(features),
    BS1 = runif(8),
    BS2 = runif(8),
    BS3 = runif(8)
  )
  ge_summarized <- summarize_by_gene_symbol(
    ge_norm,
    feature_gene_map,
    method = "mean"
  )
  expect_true(all.equal(names(ge_summarized), c("gene_symbol", "BS1", "BS2", "BS3")))
  expect_equal(length(ge_summarized$gene_symbol), length(unique(na.omit(genes))))

  expect_equal(length(setdiff(ge_summarized$gene_symbol, na.omit(genes))), 0)
  expect_equal(length(setdiff(na.omit(genes), ge_summarized$gene_symbol)), 0)
  expect_equal(
    ge_summarized[gene_symbol == "gene2"]$BS1,
    ge_norm[feature_id == "c"]$BS1
  )
  expect_equal(
    mean(ge_norm[feature_id %in% c("b", "d", "g")]$BS2),
    ge_summarized[gene_symbol == "gene1"]$BS2
  )
})

test_that("summarize_by_gene_symbol max", {
  # Genes with some duplicates to mimic features that map to multiple genes
  genes <- c(NA, "gene1", "gene2", "gene1", "gene3", "gene4", "gene1", NA, "gene2", "gene2")
  features <- letters[1:10]
  feature_gene_map <- data.table(
    featureid = features,
    genesymbol = genes
  )

  # fake gene expression data with features in a random order
  ge_norm <- data.table(
    BS1 = runif(8),
    BS2 = runif(8),
    BS3 = runif(8)
  )
  # Duplicate row 3
  ge_norm <- rbind(ge_norm, ge_norm[3, ])

  # Add another row with same mean as row3 but differen sample values
  row <- ge_norm[3, c("BS2", "BS1", "BS3")]
  names(row) <- c("BS1", "BS2", "BS3")
  ge_norm <- rbind(ge_norm, row)

  ge_norm$feature_id <- features

  ge_summarized <- summarize_by_gene_symbol(
    ge_norm,
    feature_gene_map,
    method = "max"
  )

  expect_true(all.equal(names(ge_summarized), c("gene_symbol", "BS1", "BS2", "BS3")))
  expect_equal(length(ge_summarized$gene_symbol), length(unique(na.omit(genes))))

  # duplicated genes
  expect_equal(
    ge_summarized[gene_symbol == "gene2"]$BS1,
    ge_norm[feature_id == "j"]$BS1
  )

  # Choose max
  ge_norm[, mean := rowMeans(ge_norm[, grep("^BS", names(ge_norm)), with = FALSE])]
  expect_equal(
    ge_norm[feature_id %in% c("b", "d", "g")][mean ==  max(mean)]$BS2,
    ge_summarized[gene_symbol == "gene1"]$BS2
  )

})
