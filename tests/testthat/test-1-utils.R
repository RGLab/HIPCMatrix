test_that("write_matrix", {
  output_dir <- tempdir()
  matrix_name <- "test_matrix"
  exprs <- data.table(
    feature_id = c("a", "b", "c"),
    sample1 = c(0, 0, 0),
    sample2 = c(1, 1, 1)
  )
  norm_exprs <- data.table(
    feature_id = c("a", "b", "c"),
    sample1 = c(1, 1, 1),
    sample2 = c(2, 2, 2)
  )
  sum_exprs <- data.table(
    feature_id = c("a", "b", "c"),
    sample1 = c(2, 2, 2),
    sample2 = c(1, 1, 1)
  )


  write_matrix(
    output_dir = output_dir,
    matrix_name = matrix_name,
    exprs = exprs,
    norm_exprs = norm_exprs,
    sum_exprs = sum_exprs,
    verbose = FALSE
  )

  outputs <- list.files(output_dir)
  expect_true(all(c(
    "test_matrix.tsv",
    "test_matrix.tsv.raw",
    "test_matrix.tsv.summary",
    "test_matrix.tsv.summary.orig"
  ) %in% outputs))
  expect_true(all.equal(
    fread(file.path(output_dir, "test_matrix.tsv")),
    norm_exprs
  ))
  expect_true(all.equal(
    fread(file.path(output_dir, "test_matrix.tsv.raw")),
    exprs
  ))
  expect_true(all.equal(
    fread(file.path(output_dir, "test_matrix.tsv.summary")),
    sum_exprs
  ))
  expect_true(all.equal(
    fread(file.path(output_dir, "test_matrix.tsv.summary.orig")),
    sum_exprs
  ))
})

test_that("log_message", {
  name <- "Helen"
  expect_message(
    HIPCMatrix:::log_message("Hello ", name),
    "Hello Helen"
  )
  expect_message(
    HIPCMatrix:::log_message("Hello ", name),
    "\\[\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}\\]"
  )
})

test_that("get_supp_files_dir", {
  analysis_dir <- tempdir()
  gef <- data.table(type = "PBMC", arm_accession = 1)
  supp_files_dir <- .get_supp_files_dir(analysis_dir, gef)
  expect_equal(supp_files_dir, file.path(analysis_dir, "supp_files", "PBMC_1"))

  gef <- data.table(type = "Whole blood", arm_accession = 1)
  supp_files_dir <- .get_supp_files_dir(analysis_dir, gef)
  expect_equal(supp_files_dir, file.path(analysis_dir, "supp_files", "Whole-blood_1"))
})

test_that("summarize_by_gene_symbol", {
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
    feature_gene_map
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
