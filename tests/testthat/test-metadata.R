
test_that("special case is handled correctly", {
  study <- "SDY1529"
  gef <- readRDS("test_data/sdy1529/SDY1529_day0_gef.rds")
  fas_id <- 54
  meta_data <- get_meta_data(study,
                             gef,
                             fas_id,
                             "https://www.immunespace.org")
  expect_equal(meta_data$file_location, "custom")
  expect_failure(expect_null(meta_data$custom_file_info))
  expect_equal(meta_data$custom_file_info$directory, "author_data")
  expect_equal(meta_data$custom_file_info$file_identifier_regex, "GA")
})

test_that("file location is identified correctly", {
  study <- "SDY1328"
  meta_data <- get_meta_data(study)
  expect_equal("gsm_supp_files", meta_data$file_location)

  study <- "SDY1324"
  meta_data <- get_meta_data(study)
  expect_equal("custom", meta_data$file_location)
})

test_that("custom_file_info is identified correctly", {
  study <- "SDY1324"
  meta_data <- get_meta_data(study)
  expect_failure(expect_null(meta_data$custom_file_info))
  expect_equal(meta_data$custom_file_info$directory, "raw_counts")
  expect_equal(meta_data$custom_file_info$file_identifier_regex, "RawCounts")
})
