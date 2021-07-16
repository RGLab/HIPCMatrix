mx_to_dt <- function(mx) {
  dt <- data.table(mx)
  dt$feature_id <- rownames(mx)
  dt
}

test_that("normalize_microarray works correctly", {
  affy <- readRDS("test_data/affy_raw_mx.rds")
  expect_error(
    normalize_microarray(affy),
    "max\\(exprs_mx\\) < 100"
  )
  expect_message(
    normalize_microarray(affy, force = TRUE, verbose = TRUE),
    "Forcing log2 transform"
  )
  affy_norm <- normalize_microarray(affy,
    log2_transform = FALSE
  )
  expect_equal(rownames(affy), rownames(affy_norm))
  # If quantile normalization occurred correctly each sample should have
  # the same mean
  expect_lt(diff(range(colMeans(affy_norm))), 0.001)

  illumina <- readRDS("test_data/illumina_raw_mx.rds")
  expect_error(illumina_norm <- normalize_microarray(illumina),
    regexp = NA
  )
  # Check that log2 transform was performed
  expect_lt(max(illumina_norm), 100)
})

test_that("normalize_rnaseq works correctly", {
  rna <- readRDS("test_data/rna_raw_mx.rds")
  expect_message(
    rna_norm <- normalize_rnaseq(rna, verbose = TRUE),
    "variance stabilizing transformation"
  )
  expect_equal(rownames(rna), rownames(rna_norm))
  expect_equal(colnames(rna), colnames(rna_norm))
})

test_that("normalize_matrix runs the correct mehtod", {
  library(data.table)

  affy_dt <- mx_to_dt(readRDS("test_data/affy_raw_mx.rds"))
  illumina_dt <- mx_to_dt(readRDS("test_data/illumina_raw_mx.rds"))
  rna_dt <- mx_to_dt(readRDS("test_data/rna_raw_mx.rds"))

  expect_message(
    affy_norm <- normalize_matrix(copy(affy_dt),
      "Affymetrix",
      verbose = TRUE
    ),
    "normalize_microarray"
  )
  expect_gt(max(affy_norm[, 3]), 10)
  expect_lt(max(affy_norm[, 3]), 100)
  expect_s3_class(affy_norm, "data.table")
  expect_equal(colnames(affy_norm)[1], "feature_id")

  expect_message(
    stanford_norm <- normalize_matrix(copy(affy_dt),
      "Stanford Functional Genomics Facility",
      verbose = TRUE
    ),
    "normalize_microarray"
  )
  expect_gt(max(stanford_norm[, 3]), 10)
  expect_lt(max(stanford_norm[, 3]), 100)
  expect_s3_class(stanford_norm, "data.table")
  expect_equal(colnames(stanford_norm)[1], "feature_id")

  expect_message(
    illumina_norm <- normalize_matrix(copy(illumina_dt),
      "Illumina",
      verbose = TRUE
    ),
    "log2-transform"
  )
  expect_gt(max(illumina_norm[, 3]), 10)
  expect_lt(max(illumina_norm[, 3]), 100)
  expect_s3_class(illumina_norm, "data.table")
  expect_equal(colnames(illumina_norm)[1], "feature_id")

  expect_message(
    rna_norm <- normalize_matrix(copy(rna_dt),
      "NA",
      verbose = TRUE
    ),
    "normalize_rnaseq"
  )
  expect_gt(max(rna_norm[, 3]), 10)
  expect_lt(max(rna_norm[, 3]), 100)
  expect_s3_class(rna_norm, "data.table")
  expect_equal(colnames(rna_norm)[1], "feature_id")
})

test_that("normalize_matrix handles edge cases", {
  rna <- readRDS("test_data/rna_raw_mx.rds")

  # Missing values: Should remove rows with missing values and
  # throw warning.
  mx <- rna
  mx[c(1, 5, 9), 5] <- NA
  expect_error(normalize_rnaseq(mx), "Missing values found")
  expect_error(normalize_microarray(mx), "Missing values found")
  dt <- mx_to_dt(mx)
  expect_warning(
    dt_norm <- normalize_matrix(dt, "NA"),
    "Removing 3 rows with missing values"
  )
  expect_warning(
    normalize_matrix(dt, "Affymetrix"),
    "Removing 3 rows with missing values"
  )
  expect_equal(nrow(dt) - nrow(dt_norm), 3)

  # Handles duplicated column names, returning original column names
  mx <- rna
  # duplicate first column
  mx <- cbind(mx, mx[, 1])
  colnames(mx)[ncol(mx)] <- colnames(mx)[1]
  expect_warning(mx_norm <- normalize_rnaseq(mx), "Duplicate column name: BS1000427")
  expect_true(all(colnames(mx) == colnames(mx_norm)))
  expect_warning(mx_norm <- normalize_microarray(mx), "Duplicate column name: BS1000427")
  expect_true(all(colnames(mx) == colnames(mx_norm)))
})
