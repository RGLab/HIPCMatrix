.process_illumina <- HIPCMatrix:::.process_illumina
.process_affy <- HIPCMatrix:::.process_affy
.map_feature_id_col <- HIPCMatrix:::.map_feature_id_col

test_that("map_feature_id_col correctly maps feature id column", {
  exprs_dt <- data.table(feature_id = c("A", "B"),
                         GSM1 = c(1,2),
                         GSM2 = c(3,4))
  expect_true(all.equal(exprs_dt, .map_feature_id_col(copy(exprs_dt))))
  exprs_dt <- data.table(ID_REF = c("A", "B"),
                         GSM1 = c(1,2),
                         GSM2 = c(3,4))
  expect_true(all(c("feature_id", "GSM1", "GSM2") == names(.map_feature_id_col(copy(exprs_dt)))))

  exprs_dt <- data.table(V1 = c("A", "B"),
                         GSM1 = c(1,2),
                         GSM2 = c(3,4))
  expect_true(all(c("feature_id", "GSM1", "GSM2") == names(.map_feature_id_col(copy(exprs_dt)))))

  exprs_dt <- data.table(gene_name = c("A", "B"),
                         GSM1 = c(1,2),
                         GSM2 = c(3,4))
  expect_true(all(c("feature_id", "GSM1", "GSM2") == names(.map_feature_id_col(copy(exprs_dt)))))


  exprs_dt <- data.table(rn = c("A", "B"),
                         GSM1 = c(1,2),
                         GSM2 = c(3,4))
  expect_true(all(c("feature_id", "GSM1", "GSM2") == names(.map_feature_id_col(copy(exprs_dt)))))
})

test_that("process_illumina handles subjects with no expression", {
  raw_illumina <- fread("test_data/sdy180/supp_files/Whole blood_ARM773/SDY180_raw_expression.txt")
  raw_illumina[, GSM744991.AVG_Signal.AVG_Signal := 0]
  input_file <- tempfile()
  fwrite(raw_illumina, input_file)
  processed_illumina <- .process_illumina(input_file)

  fixed_input_file <- fread(input_file)
  expect_false(any(grepl("GSM744991", names(processed_illumina))))
  expect_false(any(grepl("GSM744991", names(fixed_input_file))))
})

test_that("process_illumina handles misnamed probes", {
  raw_illumina <- fread("test_data/sdy180/supp_files/Whole blood_ARM773/SDY180_raw_expression.txt")
  raw_illumina[5, ID_REF := "BAD_PROBE_NAME"]
  input_file <- tempfile()
  fwrite(raw_illumina, input_file)
  processed_illumina <- .process_illumina(input_file)
  fixed_input_file <- fread(input_file)
  expect_false("BAD_PROBE_NAME" %in% processed_illumina[, 1])
  expect_false("BAD_PROBE_NAME" %in% fixed_input_file[, 1])
})

test_that("process_illumina returns one column per sample", {
  processed_illumina <- .process_illumina("test_data/sdy180/supp_files/Whole blood_ARM773/SDY180_raw_expression.txt")
  expect_true("GSM744994" %in% names(processed_illumina))
  expect_equal(ncol(processed_illumina), 3)
})

test_that("process_affy returns one column per sample", {
  input_files <- normalizePath(HIPCMatrix:::.select_input_files("test_data/sdy112"))
  affy_processed <- .process_affy(input_files)
  expect_equal(ncol(affy_processed), length(input_files) + 1)
  expect_equal(names(affy_processed)[1], "feature_id")
})
