test_that("normalize_microarray works correctly", {
  affy <- readRDS("test_data/affy_raw_mx.rds")
  expect_error(normalize_microarray(affy),
               "max\\(exprs_mx\\) < 100")
  expect_message(normalize_microarray(affy, force = TRUE, verbose = TRUE),
                 "Forcing log2 transform")
  affy_norm <- normalize_microarray(affy,
                                    log2_transform = FALSE)
  expect_equal(rownames(affy), rownames(affy_norm))
  # If quantile normalization occurred correctly each sample should have
  # the same mean
  expect_lt(diff(range(colMeans(affy_norm))), 0.001)

  illumina <- readRDS("test_data/illumina_raw_mx.rds")
  expect_error(illumina_norm <- normalize_microarray(illumina),
                 regexp = NA)
  # Check that log2 transform was performed
  expect_lt(max(illumina_norm), 100)
})

test_that("normalize_rnaseq works correctly", {
  rna <- readRDS("test_data/rna_raw_mx.rds")
  expect_message(rna_norm <- normalize_rnaseq(rna, verbose = TRUE),
                 "variance stabilizing transformation")
  expect_equal(rownames(rna), rownames(rna_norm))
  expect_equal(colnames(rna), colnames(rna_norm))
})

test_that("normalize_matrix works correctly", {

  library(data.table)

  affy <- readRDS("test_data/affy_raw_mx.rds")
  affy_dt <- data.table(affy)
  affy_dt$feature_id <- rownames(affy)
  illumina <- readRDS("test_data/illumina_raw_mx.rds")
  illumina_dt <- data.table(illumina)
  illumina_dt$feature_id <- rownames(illumina)
  rna <- readRDS("test_data/rna_raw_mx.rds")
  rna_dt <- data.table(rna)
  rna_dt$feature_id <- rownames(rna)

  expect_message(affy_norm <- normalize_matrix(copy(affy_dt),
                                               "Affymetrix",
                                               verbose = TRUE),
                 "normalize_microarray")
  expect_gt(max(affy_norm[, 3]), 10)
  expect_lt(max(affy_norm[, 3]), 100)
  expect_s3_class(affy_norm, "data.table")
  expect_equal(colnames(affy_norm)[1], "feature_id")

  expect_message(stanford_norm <- normalize_matrix(copy(affy_dt),
                                                   "Stanford Functional Genomics Facility",
                                                   verbose = TRUE),
                 "normalize_microarray")
  expect_gt(max(stanford_norm[, 3]), 10)
  expect_lt(max(stanford_norm[, 3]), 100)
  expect_s3_class(stanford_norm, "data.table")
  expect_equal(colnames(stanford_norm)[1], "feature_id")

  expect_message(illumina_norm <- normalize_matrix(copy(illumina_dt),
                                                   "Illumina",
                                                   verbose = TRUE),
                 "log2-transform")
  expect_gt(max(illumina_norm[, 3]), 10)
  expect_lt(max(illumina_norm[, 3]), 100)
  expect_s3_class(illumina_norm, "data.table")
  expect_equal(colnames(illumina_norm)[1], "feature_id")

  expect_message(rna_norm <- normalize_matrix(copy(rna_dt),
                                              "NA",
                                              verbose = TRUE),
                 "normalize_rnaseq")
  expect_gt(max(rna_norm[, 3]), 10)
  expect_lt(max(rna_norm[, 3]), 100)
  expect_s3_class(rna_norm, "data.table")
  expect_equal(colnames(rna_norm)[1], "feature_id")

})
