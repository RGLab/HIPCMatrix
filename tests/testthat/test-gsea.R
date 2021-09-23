test_that("gsea returns correct format", {
  eset <- SDY269$getGEMatrix("SDY269_PBMC_TIV_Geo")
  gsea_result <- gsea(eset,
    set_name = "msigdb"
  )
  expect_equal(
    colnames(gsea_result),
    c(
      "Module",
      "Coefficient",
      "NGenes",
      "Correlation",
      "Direction",
      "PValue_10log10",
      "FDR",
      "PValue"
    )
  )
  expect_setequal(
    as.character(unique(eset$study_time_collected)[2:length(unique(eset$study_time_collected))]),
    unique(gsea_result$Coefficient)
  )
})

test_that("gsea handles custom baseline", {
  eset <- SDY269$getGEMatrix("SDY269_PBMC_TIV_Geo")
  gsea_result <- gsea(eset,
    set_name = "blood_transcription",
    contrast = "study_time_collected",
    baseline = 3
  )
  expect_setequal(
    unique(gsea_result$Coefficient),
    c("7", "0")
  )
})

test_that("gsea handles custom contrast", {
  eset <- SDY269$getGEMatrix("SDY269_PBMC_TIV_Geo")
  eset$age <- sample(c("old", "young"), dim(eset)["Samples"], replace = TRUE)
  gsea_result <- gsea(eset,
    set_name = "chaussabel",
    contrast = "age"
  )
  expect_equal(unique(gsea_result$Coefficient), "young")

  gsea_result <- gsea(eset,
    set_name = "msigdb",
    contrast = "age",
    baseline = "young"
  )
  expect_equal(unique(gsea_result$Coefficient), "old")
})

test_that("run_gsea returns correct format", {
  gsea_result <- SDY269$run_gsea(matrix_name = "SDY269_PBMC_TIV_Geo")
  expect_setequal(unique(gsea_result$Coefficient), c("3 Days", "7 Days"))
  expect_lte(gsea_result$FDR[1], gsea_result$FDR[5])
  expect_lte(gsea_result$FDR[6], gsea_result$FDR[50])
})
