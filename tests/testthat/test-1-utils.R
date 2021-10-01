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
  expect_log_message(
    HIPCMatrix:::log_message("Hello ", name),
    "Hello Helen"
  )
  expect_log_message(
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
