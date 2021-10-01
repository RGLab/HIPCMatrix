### Immune Response Predictor S3 object
#' Print an ImmuneResponsePredictor object
#' @export
#' @param x an object of class ImmuneResponsePredictor
#' @param ... further arguments passed to or from other methods
print.ImmuneResponsePredictor <- function(x, ...) {
  cat("ImmuneResponsePredictor object\n")
  cat("Predictors:", nrow(x$predictors), "genes ")
  if (x$inputs$timepoint > 0) {
    cat("(fold-change at ", x$inputs$timepoint, " ", x$inputs$timepoint_unit, ")\n", sep = "")
  } else {
    cat("(expression at ", x$inputs$timepoint, " ", x$inputs$timepoint_unit, ")\n", sep = "")
  }
  cat("Response:", x$response$assay)
  cat(ifelse(x$inputs$dichotomize,
    paste0("(dichotomized: threshold = ", x$inputs$dichotomize_thresh, ")\n"),
    "\n"
  ))
  cat("Training: ", paste0(x$cohorts$training, collapse = ", "), "\n", sep = "")
  cat("Testing: ", paste0(x$cohorts$testing, collapse = ", "), "\n", sep = "")
}


# Get immune response.
#
# Returns max log fold change of \code{assay} from baseline for each
# particiapnt.
#
# See HMX$get_immune_response for full documentation
get_immune_response <- function(con,
                                assay = "hai",
                                participant_ids = NULL,
                                dichotomize = FALSE,
                                dichotomize_thresh = 4,
                                reload = FALSE) {
  args <- list(
    assay = assay,
    dichotomize = dichotomize,
    dichotomize_thresh = dichotomize_thresh
  )
  digestedArgs <- digest::digest(args)

  cache_name <- paste0("response_", digestedArgs)
  if (cache_name %in% names(con$cache) & !reload) {
    log_message("Returning response values from cache...")
    response <- con$cache[[cache_name]]
    if (!is.null(participant_ids)) {
      # Return data.table with one row per participant_id,
      # in the same order. NA values if no response data for
      # that pid
      response <- merge(data.table(participant_id = participant_ids),
        con$cache[[cache_name]],
        all.x = TRUE,
        sort = FALSE
      )
    }
    return(response)
  }


  if (!assay %in% con$availableDatasets$Name) {
    stop(sprintf("`%s` is not a valid dataset", assay))
  }


  response <- con$getDataset(assay,
    original_view = TRUE
  )

  # Get fold change from baseline
  analyte <- switch(assay,
    hai = "virus",
    neut_ab_titer = "virus",
    elisa = "analyte"
  )
  response <- response[, .(arm_accession,
    study_time_collected,
    study_time_collected_unit,
    response = value_preferred / mean(value_preferred[study_time_collected <= 0])
  ),
  by = c(analyte, "participant_id")
  ]


  # Get peak immunogenicity
  response <- response[, mean_response := mean(response, na.rm = TRUE), by = "study_time_collected"]
  response <- response[, max_response := max(mean_response), by = "arm_accession"]
  peak <- unique(response[mean_response == max_response, list(study_time_collected, arm_accession)])

  response <- merge(response, peak, by = c("study_time_collected", "arm_accession"))
  response <- response[, list(response = log2(max(response))), by = "participant_id"]
  if (dichotomize) {
    response <- response[, response := ifelse(response >= log2(dichotomize_thresh), TRUE, FALSE)]
  }
  response <- response[order(participant_id)]

  con$cache[[cache_name]] <- response


  if (!is.null(participant_ids)) {
    # Return data.table with one row per participant_id,
    # in the same order. NA values if no response data for
    # that pid
    response <- merge(data.table(participant_id = participant_ids),
      con$cache[[cache_name]],
      all.x = TRUE,
      sort = FALSE
    )
  }

  return(response)
}



# Get genes differentially expressed over time
# See HMX$get_de_genes
get_de_genes <- function(con,
                         timepoint,
                         fc_thresh = 0.58,
                         timepoint_unit = "Days",
                         cohorts = NULL,
                         reload = FALSE) {

  # Check cache
  args <- list(
    timepoint = timepoint,
    fc_thresh = fc_thresh,
    timpeoint_unit = timepoint_unit,
    cohorts = cohorts
  )
  digestedArgs <- digest::digest(args)

  if (!"de_genes" %in% names(con$cache)) con$cache$de_genes <- list()
  if (digestedArgs %in% names(con$cache$de_genes) & !reload) {
    log_message("Returning differentially expressed genes from cache...")
    return(con$cache$de_genes[[digestedArgs]]$data)
  }


  # Validate params

  if (!is.null(cohorts)) {
    bad_cohorts <- !cohorts %in% con$listGEMatrices()$cohort_type
    if (any(bad_cohorts)) {
      stop(
        "'", cohorts[bad_cohorts], "' is not a valid cohort. \n",
        "Choose from: ", paste0(con$listGEMatrices()$cohort_type,
          collapse = ", "
        )
      )
    }
  }

  if (!tolower(paste(timepoint, timepoint_unit)) %in% tolower(con$listGEAnalysis()$coefficient)) {
    stop(
      timepoint, " ", timepoint_unit, " is not a valid timepoint. Please choose from ",
      paste0(con$listGEAnalysis()$coefficient, collapse = ", ")
    )
  }


  colfilter <- makeFilter(c(
    "coefficient",
    "EQUAL",
    paste(timepoint, timepoint_unit)
  ))

  if (!is.null(cohorts)) {
    cohorts_filter <- makeFilter(c(
      "cohort",
      "IN",
      paste(cohorts,
        collapse = ";"
      )
    ))
    colfilter <- rbind(cohorts_filter, colfilter)
  }

  dgear <- con$getGEAnalysis(colFilter = colfilter)
  setnames(dgear, "log_fold-change", "log_fc")
  dgear <- dgear[abs(log_fc) > fc_thresh]
  de_genes <- unique(dgear$gene_symbol)


  con$cache$de_genes[[digestedArgs]] <- list(
    data = de_genes,
    args = args
  )

  de_genes
}


get_fc_mx <- function(eset,
                      timepoint,
                      timepoint_unit = "Days",
                      features = NULL) {

  # Get only subject with requested timepoint OR baseline AND appropriate unit
  eset <- eset[, eset$study_time_collected %in% c(0, timepoint) &
    tolower(eset$study_time_collected_unit) == tolower(timepoint_unit)]
  pd <- data.table(Biobase::pData(eset))

  # Calculate fold-change at selected timepoint
  #
  # In some cases, e.g. SDY180, there may be multiple observations per
  # participant_id + timePt, therefore need to arbitrarily select the first
  # observation prior to ensuring subject has BOTH requested timePt and
  # baseline.
  setkeyv(pd, c("study_time_collected", "participant_id"))
  pd <- pd[, .SD[1L], by = key(pd)]

  if (timepoint > 0) {
    pids <- pd[, .N, by = participant_id][N > 1, participant_id]
    pd <- pd[participant_id %in% pids]
    pd <- pd[order(participant_id, study_time_collected)]
    eset <- eset[, pd$biosample_accession]

    # Calculate fold-change
    later_tp <- Biobase::exprs(eset[, pd[study_time_collected == timepoint, biosample_accession]])
    baseline <- Biobase::exprs(eset[, pd[study_time_collected == 0, biosample_accession]])
    FC <- t(later_tp - baseline)
  } else {
    pids <- pd$participant_id
    pd <- pd[order(participant_id, study_time_collected)]
    eset <- eset[, pd$biosample_accession]
    FC <- t(Biobase::exprs(eset))
  }

  rownames(FC) <- pd[match(rownames(FC), pd$biosample_accession), participant_id]

  if (!is.null(features)) {
    # Add NA cols for features not present
    if (any(!features %in% colnames(FC))) {
      new_features <- features[!features %in% colnames(FC)]
      FC <- cbind(
        FC,
        matrix(
          nrow = nrow(FC),
          ncol = length(new_features),
          dimnames = list(rownames(FC), new_features)
        )
      )
    }
    FC <- FC[, features]
  }

  FC
}


# Given a matrix of fold-change for each feature for each observation,
# and a vector of responses, return a vector of genes which are determined
# to be predictors using elastic net pentalty.
select_features <- function(FC,
                            response_vector,
                            dichotomize = FALSE) {

  # Check prerequisites for glmnet
  if (is.null(dim(FC))) {
    stop("Training cohort results in one participant at peak time point.
         This is insufficient for performing glmnet function.
         Please select different training cohort and retry.")
  }

  if (ncol(FC) <= 1) {
    stop("At least two genes are required for prediction.
             Try lowering or diabling the log fold-change filter in the
             'Additional options' section to increase the number of selected
             genes.")
  }

  if (length(unique(response_vector)) == 1) {
    stop("All response calls in training cohort are the same,
         which causes the glmnet() to fail as 'y' is constant and the null
         deviance cannot be calculated in order to standardize the data.
         Please select a different training cohort.")
  }

  if (nrow(FC) != length(response_vector)) {
    stop(
      "Number of observations in response_vector (",
      length(response_vector),
      ") not equal to number of rows in FC (",
      nrow(FC), ")"
    )
  }

  if (isTRUE(dichotomize) & class(response_vector) != "logical") {
    stop("response_vector must be logical when dichotomize = TRUE")
  }

  # Elastic net

  if (dichotomize) {
    # alpha = 0.5 means equal ridge and lasso penalty terms
    fit <- glmnet::glmnet(FC, as.factor(response_vector), alpha = 0.5, family = "binomial")
  } else {
    fit <- glmnet::glmnet(FC, response_vector, alpha = 0.5)
  }

  cv_fit <- glmnet::cv.glmnet(FC, response_vector)
  # use lambda.min to get lambda that resulted in smallest sum
  coef <- stats::predict(fit, s = cv_fit$lambda.min, type = "coefficients")


  selected_features <- names(which(abs(coef[, 1]) > 0))
  selected_features <- grep("Intercept", selected_features, invert = TRUE, value = TRUE)
  if (length(selected_features) < 2) {
    opts_chunk$set(eval = FALSE, cache = FALSE)
    stop("No features were selected as predictive. You may try to remove the fold change filtering under 'Additional options'.")
  }

  # Need more obs than features
  nFeatures <- length(selected_features)
  nObs <- length(response_vector)
  if (nObs <= nFeatures) {
    log_message(
      "You selected as many or more features (",
      nFeatures,
      ") than observations (",
      nObs,
      ").\nThe",
      nObs - 2,
      "most significant features will be kept.\n"
    )
    selected_features <- names(sort(coef[selected_features, ])[1:(nObs - 2)])
  }

  selected_features
}


# Fit linear model with selected_features as predictors
get_fit <- function(FC,
                    response_vector,
                    selected_features,
                    dichotomize = FALSE) {

  # Lasso: Get coefficients
  FC <- FC[, selected_features]
  FC <- FC[, selected_features]

  # lasso #IN FC, dichotomize, selected_features #OUT predictor_table
  FC <- data.frame(FC, check.names = FALSE)
  formula <- stats::as.formula(paste0("outcome~`", paste(colnames(FC), collapse = "`+`"), sep = "`"))
  FC$outcome <- response_vector
  if (dichotomize) {
    fit <- stats::glm(formula, FC, family = "binomial")
  } else {
    fit <- stats::lm(formula, FC)
  }

  fit
}


# Train immune response predictors
#
# See HMX$train_immune_response_predictors for full documentation

train_immune_response_predictors <- function(con,
                                             cohorts,
                                             timepoint,
                                             assay = "hai",
                                             timepoint_unit = "Days",
                                             use_only_de_genes = timepoint > 0,
                                             fc_thresh = 0.58,
                                             dichotomize = FALSE,
                                             dichotomize_thresh = 4,
                                             reload = FALSE) {
  # Check cache
  args <- list(
    cohorts = cohorts,
    timepoint = timepoint,
    assay = assay,
    timepoint_unit = timepoint_unit,
    use_only_de_genes = use_only_de_genes,
    fc_thresh = fc_thresh,
    dichotomize = dichotomize,
    dichotomize_thresh = dichotomize_thresh
  )

  digestedArgs <- digest::digest(args)
  hashes <- sapply(con$immune_response_predictors, `[[`, "hash")
  if (digestedArgs %in% hashes & !reload) {
    log_message("returning model fit from cache.")
    return(which(hashes == digestedArgs))
  }


  # Validate params

  bad_cohorts <- !cohorts %in% con$listGEMatrices()$cohort_type
  if (any(bad_cohorts)) {
    stop(
      "'", cohorts[bad_cohorts], "' is not a valid cohort. \n",
      "Choose from: ", paste0(con$listGEMatrices()$cohort_type,
        collapse = ", "
      )
    )
  }

  if (!tolower(paste(timepoint, timepoint_unit)) %in% tolower(con$listGEAnalysis()$coefficient) &
    timepoint > 0) {
    stop(
      timepoint, " ", timepoint_unit, " is not a valid timepoint. Please choose from ",
      paste0(unique(con$listGEAnalysis()$coefficient), collapse = ", ")
    )
  }
  eset <- con$getGEMatrix(cohortType = cohorts)

  # Subset to DE genes
  if (use_only_de_genes) {
    log_message("Subsetting to differentially expressed genes...")
    de_genes <- con$get_de_genes(
      timepoint = timepoint,
      fc_thresh = fc_thresh,
      timepoint_unit = timepoint_unit,
      cohorts = cohorts
    )
    log_message(length(de_genes), " differentially expressed genes found.")
    eset <- eset[Biobase::featureNames(eset) %in% de_genes, ]
  }

  # Get matrix of fold-change from day 0 to timepoint.
  # If timepoint = 0, then just return expression values.
  log_message(
    "Deriving fold-change from baseline to ",
    timepoint,
    " ",
    timepoint_unit,
    "..."
  )
  FC <- get_fc_mx(
    eset,
    timepoint,
    timepoint_unit
  )

  # Get immune response for each participant.
  log_message("Calculating immune response for ", assay, "...")
  response <- con$get_immune_response(
    assay = assay,
    participant_ids = rownames(FC),
    dichotomize = dichotomize,
    dichotomize_thresh = dichotomize_thresh
  )

  # remove any participants with NA values
  response <- response[!is.na(response)]
  FC <- FC[rownames(FC) %in% unique(response$participant_id), ]
  response_vector <- response$response

  if (timepoint > 0) {
    if (!use_only_de_genes & fc_thresh > 0) { # i.e: not using the GEAR but using a FC thresh on the matrix
      FC <- FC[, log(apply(abs(FC), 2, max)) > fc_thresh, drop = FALSE]
    }
  }

  # Select features using elastic net
  log_message("Selecting features using elastic net...")
  features <- select_features(
    FC = FC,
    response_vector = response_vector,
    dichotomize = dichotomize
  )
  log_message(length(features), " features selected.")


  # Fit linear model using selected features
  log_message("Finding model fit using selected features...")
  fit <- get_fit(FC,
    response_vector,
    features,
    dichotomize = dichotomize
  )

  irp <- list(
    hash = digestedArgs,
    model = fit,
    cohorts = list(
      training = cohorts,
      testing = NULL
    ),
    predictors = format_predictor_table(fit),
    response = list(
      assay = assay,
      dichotomize = dichotomize,
      dichotomize_thresh = dichotomize_thresh,
      data = list(
        training = response[, .(
          participant_id,
          observed = response,
          predicted = fit$fitted.values[participant_id],
          cohort = eset[, match(participant_id, eset$participant_id)]$cohort_type
        )],
        testing = NULL
      )
    ),
    inputs = args,
    FC = FC[, features]
  )

  class(irp) <- "ImmuneResponsePredictor"

  irp_index <- length(con$immune_response_predictors) + 1
  con$immune_response_predictors[[irp_index]] <- irp

  invisible(irp_index)
}

# by default, run on all irp objects
predict_response <- function(con,
                             cohort,
                             irp_index) {

  # TODO: Predict differently when dichotomize = TRUE?

  if (is.null(irp_index) | is.null(con$immune_response_predictors[[irp_index]])) {
    stop("Invalid irp_index")
  }
  # Validate params
  if (length(cohort) != 1) stop("one cohort must be specified.")
  bad_cohorts <- !cohort %in% con$listGEMatrices()$cohort_type
  if (any(bad_cohorts)) {
    stop(
      "'", cohort, "' is not a valid cohort. \n",
      "Choose from: ", paste0(con$listGEMatrices()$cohort_type,
        collapse = ", "
      )
    )
  }

  irp <- con$get_irp(irp_index)

  timepoint <- irp$inputs$timepoint
  timepoint_unit <- irp$inputs$timepoint_unit

  eset <- con$getGEMatrix(cohortType = cohort)

  # Get matrix of fold-change from day 0 to timepoint.
  # If timepoint = 0, then just return expression values.
  log_message(
    "Deriving fold-change from baseline to ",
    timepoint,
    " ",
    timepoint_unit,
    "..."
  )

  FC <- get_fc_mx(
    eset,
    timepoint,
    timepoint_unit,
    features = irp$predictors$gene_symbol
  )

  log_message("Deriving predicted values...")
  newdata <- data.frame(FC)
  colnames(newdata) <- colnames(FC)
  rownames(newdata) <- rownames(FC)
  predicted_values <- stats::predict(irp$model,
    newdata = newdata,
    type = "response",
    na.action = stats::na.exclude
  )

  new_FC <- rbind(
    con$get_irp(irp_index)$FC,
    FC[, colnames(con$get_irp(irp_index)$FC)]
  )
  new_FC <- new_FC[!duplicated(new_FC), ]
  con$immune_response_predictors[[irp_index]]$FC <- new_FC

  predicted_values
}


# Test existing irp model on a different cohort
# see HMS$test_immune_response_predictors for full documentation
test_immune_response_predictors <- function(con,
                                            cohorts,
                                            irp_index = NULL,
                                            reload = FALSE) {
  if (length(con$immune_response_predictors) == 0) {
    stop("No irp objects found in con. First run con$train_immune_response_predictors()")
  }

  if (is.null(irp_index)) {
    log_message("irp_index not specified. Using most recent...")
    irp_index <- length(con$immune_response_predictors)
  }

  if (irp_index > length(con$immune_response_predictors)) {
    stop("Invalid irp_index")
  }

  irp <- con$get_irp(irp_index)

  assay <- irp$inputs$assay
  dichotomize <- irp$inputs$dichotomize
  dichotomize_thresh <- irp$inputs$dichotomize_thresh

  res <- lapply(cohorts, function(cohort_test) {

    # When reload = TRUE, run all cohorts not in TRAINING set.
    # when reload = FALSE, only run cohorts not already in TESTING set.
    if (cohort_test %in% irp$cohorts$training) {
      result <- irp$response$data$training[cohort == cohort_test]
    } else if (cohort_test %in% irp$cohorts$testing & !reload) {
      result <- irp$response$data$testing[cohort == cohort_test]
    } else {
      predicted_values <- con$predict_response(
        cohort = cohort_test,
        irp_index = irp_index
      )

      result <- con$get_immune_response(
        assay = assay,
        participant_ids = names(predicted_values),
        dichotomize = dichotomize,
        dichotomize_thresh = dichotomize_thresh
      )
      setnames(result, "response", "observed")
      result[, predicted := predicted_values[participant_id]]
      result[, cohort := cohort_test]

      # Don't add to testing list if already there ie reload = TRUE
      if (!cohort_test %in% irp$cohorts$testing) {
        con$immune_response_predictors[[irp_index]]$response$data$testing <- rbind(
          irp$response$data$testing,
          result
        )
        con$immune_response_predictors[[irp_index]]$cohorts$testing <- c(irp$cohorts$testing, cohort_test)
      } else {
        # rm previous entry and add new
        con$immune_response_predictors[[irp_index]]$response$data$testing <- irp$response$data$testing[cohort != cohort_test]
        con$immune_response_predictors[[irp_index]]$response$data$testing <- rbind(
          irp$response$data$testing,
          result
        )
      }
    }

    result[, set := ifelse(cohort_test %in% irp$cohorts$training,
      "Training",
      "Testing"
    )]

    result
  })

  rbindlist(res)
}

# Run IRP
run_irp <- function(con,
                    cohorts_train,
                    cohorts_test,
                    timepoint,
                    use_only_de_genes,
                    assay,
                    timepoint_unit,
                    fc_thresh,
                    dichotomize,
                    dichotomize_thresh,
                    reload) {
  irp_index <- con$train_immune_response_predictors(
    cohorts = cohorts_train,
    timepoint = timepoint,
    assay = assay,
    timepoint_unit = timepoint_unit,
    use_only_de_genes = use_only_de_genes,
    fc_thresh = fc_thresh,
    dichotomize = dichotomize,
    dichotomize_thresh = dichotomize_thresh,
    reload = reload
  )

  con$test_immune_response_predictors(
    cohorts = cohorts_test,
    irp_index = irp_index,
    reload = reload
  )

  con$get_irp(irp_index)
}

#' Format Predictor Table
#'
#' Create nicely formatted table of predictors from fit result.
#'
#' @param fit fit object from \code{glm} or \code{lm}
#'
#' @export
format_predictor_table <- function(fit) {
  sum_fit <- summary(fit)
  sum_fit_coef <- sum_fit$coefficients
  # Strip beginning and ending "`" which get added to names with weird symbols
  rownames(sum_fit_coef) <- gsub("^`|`$", "", rownames(sum_fit_coef))
  pred_cIdx <- grep("value|Pr", colnames(sum_fit_coef))
  predictor_table <- sum_fit_coef[, pred_cIdx][-1, ]
  colnames(predictor_table) <- c("statistic", "p-value")


  predictor_table <- data.table(cbind(
    data.table(selected_features = rownames(predictor_table)),
    predictor_table
  ))
  predictor_table[, .(
    gene_symbol = selected_features,
    statistic = as.numeric(statistic),
    `p-value` = as.numeric(`p-value`)
  )]
}
