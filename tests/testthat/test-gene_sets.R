test_that("update_geneset_symbols returns correct format", {
  btm_updated <- update_geneset_symbols("orig_btm_list")
  expect_equal(names(btm_updated), c("pathway", "SYMBOL"))
  expect_s3_class(btm_updated, "data.table")

  emory_updated <- update_geneset_symbols("orig_emory")
  expect_equal(names(emory_updated), c("pathway", "SYMBOL"))
  expect_s3_class(emory_updated, "data.table")
})

test_that("update_geneset_symbols result matches package data", {
  btm_updated <- update_geneset_symbols("orig_btm_list")
  expect_true(all.equal(btm_updated, HIPCMatrix::updated_btm_df))
})
