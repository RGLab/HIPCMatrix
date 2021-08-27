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

test_that("HMX$runGEAnalysis returns correct object", {
  de_result <- con$runGEAnalysis()
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
})

# test_that("HMX$uploadGEAnalysisResults", {
#   with_mock_api({
#     con$uploadGEAnalysisResults()
#   })
# })
