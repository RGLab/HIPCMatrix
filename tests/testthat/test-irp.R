test_that("get_immune_response handles bad inputs", {
  expect_error(
    SDY269$get_immune_response(assay = "assay"),
    "`assay` is not a valid dataset"
  )
})

test_that("get_immune_response returns correct format", {
  resp <- SDY269$get_immune_response()
  expect_equal(colnames(resp), c("participant_id", "response"))
  expect_type(resp$response, "double")
  expect_equal(length(unique(resp$participant_id)), nrow(resp))
})

test_that("get_immune_response with dichotomize returns correct format", {
  resp_dich <- SDY269$get_immune_response(dichotomize = TRUE)
  expect_equal(colnames(resp_dich), c("participant_id", "response"))
  expect_type(resp_dich$response, "logical")
  expect_equal(length(unique(resp_dich$participant_id)), nrow(resp_dich))
})

test_that("get_immune_response only returns correct pids", {
  demog <- SDY269$getDataset("demographics")
  pids <- demog$participant_id[1:10]
  resp <- SDY269$get_immune_response(participant_ids = pids)
  expect_equal(pids, resp$participant_id)

  # pid which is not present in assay should still show up, with NA value
  pids <- c(pids, "TEST")
  resp <- SDY269$get_immune_response(participant_ids = pids)
  expect_equal(pids, resp$participant_id)
  expect_true(is.na(resp[participant_id == "TEST", response]))
})

test_that("get_immune_response neut_ab_titer", {
  con <- HMX$new("SDY180")
  resp <- con$get_immune_response("neut_ab_titer")
  expect_equal(colnames(resp), c("participant_id", "response"))
  expect_type(resp$response, "double")
  expect_equal(length(unique(resp$participant_id)), nrow(resp))
})

test_that("get_immune_response elisa", {
  resp <- SDY269$get_immune_response(assay = "elisa")
  expect_equal(colnames(resp), c("participant_id", "response"))
  expect_type(resp$response, "double")
  expect_equal(length(unique(resp$participant_id)), nrow(resp))
})

test_that("get_de_genes filters correctly", {
  de_genes <- SDY269$get_de_genes(timepoint = 7)

  de_genes_tiv <- SDY269$get_de_genes(timepoint = 7, cohorts = "TIV Group 2008_PBMC")
  de_genes_laiv <- SDY269$get_de_genes(timepoint = 7, cohorts = "LAIV group 2008_PBMC")

  expect_gt(length(de_genes), length(de_genes_tiv))
  expect_equal(de_genes, union(de_genes_tiv, de_genes_laiv))
})

test_that("get_de_genes handles bad inputs", {
  expect_error(
    SDY269$get_de_genes(SDY269, cohorts = "fake cohort"),
    "'fake cohort' is not a valid cohort."
  )

  expect_error(
    SDY269$get_de_genes(SDY269, timepoint = 9),
    "9 Days is not a valid timepoint"
  )
})


select_features <- HIPCMatrix:::select_features
test_that("select_features handles bad inputs", {
  # only one subject at selected timepoint
  FC <- rnorm(30)
  response_vector <- sample(0:6, 1)
  expect_error(
    select_features(FC, response_vector),
    "Training cohort results in one participant at peak time point"
  )

  # Only one gene selected
  FC <- matrix(rnorm(30), nrow = 30)
  response_vector <- sample(0:6, 30, replace = TRUE)
  expect_error(
    select_features(FC, response_vector),
    "At least two genes are required for prediction"
  )

  # All same responses
  FC <- matrix(rnorm(30 * 50), nrow = 30, ncol = 50)
  response_vector <- rep(1, 30)
  expect_error(
    select_features(FC, response_vector),
    "All response calls in training cohort are the same"
  )

  # Mismatched dimensions
  FC <- matrix(rnorm(30 * 50), nrow = 30, ncol = 50)
  response_vector <- sample(0:6, 35, replace = TRUE)
  expect_error(
    select_features(FC, response_vector),
    "Number of observations in response_vector \\(35\\) not equal to number of rows in FC \\(30\\)"
  )

  # data type for dichotomize
  FC <- matrix(rnorm(30 * 50), nrow = 30, ncol = 50)
  response_vector <- sample(0:6, 30, replace = TRUE)
  expect_error(
    select_features(FC, response_vector, dichotomize = TRUE),
    "response_vector must be logical when dichotomize = TRUE"
  )
})

test_that("select_features returns correct format", {
  FC <- matrix(rnorm(30 * 50), nrow = 30, ncol = 50)
  rownames(FC) <- paste0("sample", 1:30)
  colnames(FC) <- paste0("gene", 1:50)
  # Make response dependent on genes 1:10
  response_vector <- rowSums(FC[, 1:10])
  features <- select_features(FC, response_vector)
  expect_type(features, "character")
  expect_true(all(grepl("gene", features)))

  # Same but dichotomize
  response_vector <- response_vector >= 0
  features <- select_features(FC, response_vector)
  expect_type(features, "character")
  expect_true(all(grepl("gene", features)))
})

test_that("select_features selects fewer features than obs", {
  set.seed(42)
  FC <- matrix(rnorm(5 * 50), nrow = 5, ncol = 50)
  rownames(FC) <- paste0("sample", 1:5)
  colnames(FC) <- paste0("gene", 1:50)
  # Make response dependent on genes 1:10
  response_vector <- rowSums(FC[, 1:10])
  expect_message(
    features <- select_features(FC, response_vector),
    "You selected as many or more features"
  )
  expect_length(features, 3)
})

get_fit <- HIPCMatrix:::get_fit
test_that("get_fit returns correct format", {
  FC <- matrix(rnorm(30 * 50), nrow = 30, ncol = 50)
  rownames(FC) <- paste0("sample", 1:30)
  colnames(FC) <- paste0("gene", 1:50)
  # Make response dependent on genes 1:10
  response_vector <- rowSums(FC[, 1:10])
  response_vector <- response_vector + rnorm(30)
  features <- colnames(FC)[1:10]
  fit <- get_fit(FC, response_vector, features)
  expect_s3_class(fit, "lm")
  expect_setequal(names(fit$coefficients), c("(Intercept)", features))
  # coefs should all be close to 1
  expect_gt(mean(fit$coefficients[2:length(fit$coefficients)]), 0.9)
  expect_lt(mean(fit$coefficients[2:length(fit$coefficients)]), 1.1)

  # dichotomize
  response_vector <- response_vector > 0
  # add a bit more randomness so it does not fit perfectly...
  response_sample <- sample(length(response_vector), 5)
  response_vector[response_sample] <- !response_vector[response_sample]
  fit <- get_fit(FC, response_vector, features, dichotomize = TRUE)
  expect_s3_class(fit, "glm")
  expect_setequal(names(fit$coefficients), c("(Intercept)", features))
})


test_that("get_immune_response_predictor works", {
  pred7 <- SDY269$train_immune_response_predictors(
    cohorts = "LAIV group 2008_PBMC",
    timepoint = 7
  )

  pred_resp <- SDY269$predict_response(
    cohorts = "TIV Group 2008_PBMC",
    timepoint = 7,
    fit = pred7
  )
  # Should now have de_genes and response in cache
  expect_true("de_genes" %in% names(SDY269$cache))
  expect_length(grep("irp_fit_", names(SDY269$cache)), 1)

  expect_s3_class(pred7, "lm")
})

test_that("test_immune_response_predictors returns correct format", {
  fit <- SDY269$train_immune_response_predictors(
    cohorts = "LAIV group 2008_PBMC",
    timepoint = 7
  )

  result_train <- SDY269$test_immune_response_predictors(
    "LAIV group 2008_PBMC",
    timepoint = 7,
    fit = fit
  )

  result_test <- SDY269$test_immune_response_predictors(
    "TIV Group 2008_PBMC",
    timepoint = 7,
    fit = fit
  )

  expect_equal(colnames(result_train), c("participant_id", "observed", "predicted"))
  expect_equal(colnames(result_test), c("participant_id", "observed", "predicted"))
  preds <- result_train$predicted
  names(preds) <- result_train$participant_id
  expect_mapequal(preds, fit$fitted.values)
})

#
# test_that("run IRP", {
#
#   cohorts_train <- "LAIV group 2008_PBMC"
#   cohorts_test <- "TIV Group 2008_PBMC"
#   timepoint = 7
#   assay <- "elisa"
#
#   fit <- SDY269$train_immune_response_predictors(cohorts_train,
#                                                  timepoint,
#                                                  assay = assay)
#
#   result_train <- SDY269$test_immune_response_predictors(cohorts_train,
#                                                          7,
#                                                          fit,
#                                                          assay = assay)
#   results_test <- SDY269$test_immune_response_predictors(cohorts_test,
#                                                          7,
#                                                          fit,
#                                                          assay = assay)
#
# })
