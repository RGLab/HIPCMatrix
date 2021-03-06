eset <- SDY269$getGEMatrix("SDY269_PBMC_TIV_Geo",
  outputType = "normalized",
  annotation = "latest"
)

test_that("find_de_genes_eBayes returns correct format", {
  de_result <- find_de_genes_eBayes(eset)

  # List with results for each non-0 timepoint
  expect_equal(names(de_result), c("coef3 Days", "coef7 Days"))
  expect_s3_class(de_result[[1]], "data.table")
  expect_equal(names(de_result[[1]]), names(de_result[[2]]))
  expect_equal(
    names(de_result[[1]]),
    c(
      "FeatureId",
      "gene_symbol",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B",
      "coefficient"
    )
  )
})

test_that("find_de_genes_eBayes uses correct p-val cutoff", {
  de_result <- find_de_genes_eBayes(eset, p_val_cutoff = 0.02)
  expect_lt(max(de_result[[1]]$adj.P.Val), 0.02)


  de_result <- find_de_genes_eBayes(eset, p_val_cutoff = 0.05)
  expect_lt(max(de_result[[1]]$adj.P.Val), 0.05)

  # Use lowest 100 if < 100 genes under cutoff
  de_result <- find_de_genes_eBayes(eset, p_val_cutoff = 0.001)
  expect_equal(nrow(de_result[[1]]), 100)
})

# test_that("find_de_genes_deseq works", {
#   con <- HMX$new("SDY1256")
#   eset_norm <- con$getGEMatrix("SDY1256_WholeBlood_EPIC001_geo")
#   expect_error(
#     de_result <- find_de_genes_deseq(eset),
#     "some values in assay are not integers"
#   )
#   rm(eset_norm)
# })

test_that("HMX$run_de_analysis returns correct object", {
  de_result <- SDY269$run_de_analysis()
  expect_s3_class(de_result, "data.table")
  expect_equal(
    names(de_result),
    c(
      "feature_id",
      "gene_symbol",
      "log_fc",
      "ave_expr",
      "statistic",
      "p_value",
      "adj_p_val",
      "B",
      "analysis_accession",
      "coefficient"
    )
  )

  sdy28 <- HMX$new("SDY28")
  expect_log_message(
    de_result <- sdy28$run_de_analysis(),
    "Insufficient data"
  )
  expect_true(is.null(de_result))

  sdy406 <- HMX$new("SDY406")
  expect_log_message(
    de_result <- sdy406$run_de_analysis(),
    "Insufficient data"
  )
})


test_that("get_de_compatible_runs returns correct format", {
  implied_de <- SDY269$get_de_compatible_runs()
  expect_equal(nrow(implied_de), 4)
})
# test_that("HMX$uploadGEAnalysisResults", {
#   with_mock_api({
#     con$uploadGEAnalysisResults()
#   })
# })
