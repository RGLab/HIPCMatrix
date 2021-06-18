test_that("write_matrix", {
  output_dir <- tempdir()
  matrix_name <- "test_matrix"
  exprs <- data.table(feature_id = c("a", "b", "c"),
                      sample1 = c(0, 0, 0),
                      sample2 = c(1, 1, 1))
  norm_exprs <- data.table(feature_id = c("a", "b", "c"),
                           sample1 = c(1,1,1),
                           sample2 = c(2,2,2))
  sum_exprs <- data.table(feature_id = c("a", "b", "c"),
                          sample1 = c(2,2,2),
                          sample2 = c(1,1,1))


  write_matrix(output_dir = output_dir,
               matrix_name = matrix_name,
               exprs = exprs,
               norm_exprs = norm_exprs,
               sum_exprs = sum_exprs,
               verbose = FALSE)

  outputs <- list.files(output_dir)
  expect_true(all(c(
    "test_matrix.tsv",
    "test_matrix.tsv.raw",
    "test_matrix.tsv.summary",
    "test_matrix.tsv.summary.orig"
  ) %in% outputs))
  expect_true(all.equal(fread(file.path(output_dir, "test_matrix.tsv")),
                        norm_exprs))
  expect_true(all.equal(fread(file.path(output_dir, "test_matrix.tsv.raw")),
                        exprs))
  expect_true(all.equal(fread(file.path(output_dir, "test_matrix.tsv.summary")),
                        sum_exprs))
  expect_true(all.equal(fread(file.path(output_dir, "test_matrix.tsv.summary.orig")),
                        sum_exprs))
})


