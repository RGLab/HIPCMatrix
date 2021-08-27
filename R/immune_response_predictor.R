
#' Get immune response.
#'
#' Returns max log fold change of \code{assay} from baseline for each
#' particiapnt.
#'
#' @param con HMX connection
#' @param assay "hai", "neut_ab_response", or "elisa"
#' @param participant_ids character vector of participant_ids
#' @param dichotomize Dichotomize result? If FALSE, max log fold change is
#' returned. If TRUE, returns TRUE if max log fold change is greater than
#' \code{dichotomize_value}
#' @param dichotomize_thresh Value to use for dichotomizing result if
#' \code{dichotomize} is TRUE.
#'
get_immune_response <- function(con,
                                assay = "hai",
                                participant_ids = NULL,
                                dichotomize = FALSE,
                                dichotomize_thresh = 4) {
  if (!assay %in% con$availableDatasets$Name) {
    stop(sprintf("`%s` is not a valid dataset", assay))
  }


  response <- con$getDataset(assay,
    original_view = TRUE
  )
  if (!is.null(participant_ids)) {
    response <- response[participant_id %in% participant_ids]
  }

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
  return(response)
}



#' Get genes differentially expressed over time
#'
#' @param con HMX connection
#' @param timepoint Timepoint for finding differentially expressed genes
#' @param fc_thresh fold-change threshold for determining whether genes
#' are differentially expressed
#' @param timepoint_unit "Days" or "Hours"
#' @param cohorts character vector of cohorts to include
#'
get_de_genes <- function(con,
                         timepoint,
                         fc_thresh = 0.58,
                         timepoint_unit = "Days",
                         cohorts = NULL) {


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

  de_genes
}


get_fc_mx <- function(eset,
                      timepoint,
                      timepoint_unit = "Days") {

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
  if (timepoint > 0) {
    setkeyv(pd, c("study_time_collected", "participant_id"))
    pd <- pd[, .SD[1L], by = key(pd)]
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
    message(
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
  formula <- as.formula(paste0("outcome~`", paste(colnames(FC), collapse = "`+`"), sep = "`"))
  FC$outcome <- response_vector
  if (dichotomize) {
    relasso <- glm(formula, FC, family = "binomial")
  } else {
    relasso <- lm(formula, FC)
  }

  relasso
}

#' Get Immune Response Predictors
#'
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
  cache_name <- paste0("irp_fit_", digestedArgs)
  if (cache_name %in% names(con$cache) & !reload) {
    return(con$cache[[cache_name]])
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

  if (!tolower(paste(timepoint, timepoint_unit)) %in% tolower(con$listGEAnalysis()$coefficient)) {
    stop(
      timepoint, " ", timepoint_unit, " is not a valid timepoint. Please choose from ",
      paste0(con$listGEAnalysis()$coefficient, collapse = ", ")
    )
  }
  eset <- con$getGEMatrix(cohortType = cohorts)

  # Subset to DE genes
  if (use_only_de_genes) {
    de_genes <- con$get_de_genes(
      timepoint = timepoint,
      fc_thresh = fc_thresh,
      timepoint_unit = timepoint_unit,
      cohorts = cohorts
    )
    eset <- eset[Biobase::featureNames(eset) %in% de_genes, ]
  }

  # Get matrix of fold-change from day 0 to timepoint.
  # If timepoint = 0, then just return expression values.
  FC <- get_fc_mx(
    eset,
    timepoint,
    timepoint_unit
  )

  # Get immune response for each participant.
  response <- con$get_immune_response(
    assay = assay,
    participant_ids = rownames(FC),
    dichotomize = dichotomize
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
  features <- select_features(
    FC = FC,
    response_vector = response_vector,
    dichotomize = dichotomize
  )

  # Fit linear model using selected features
  fit <- get_fit(FC,
    response_vector,
    features,
    dichotomize = dichotomize
  )

  con$cache[[cache_name]] <- fit

  fit
}

test_immune_response_predictors <- function(con,
                                            cohorts,
                                            timepoint,
                                            fit,
                                            assay = "hai",
                                            timepoint_unit = "Days",
                                            fc_thresh = 0.58,
                                            dichotomize = FALSE,
                                            dichotomize_thresh = 4) {
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

  if (!tolower(paste(timepoint, timepoint_unit)) %in% tolower(con$listGEAnalysis()$coefficient)) {
    stop(
      timepoint, " ", timepoint_unit, " is not a valid timepoint. Please choose from ",
      paste0(con$listGEAnalysis()$coefficient, collapse = ", ")
    )
  }
  eset <- con$getGEMatrix(cohortType = cohorts)

  # Get matrix of fold-change from day 0 to timepoint.
  # If timepoint = 0, then just return expression values.
  FC <- get_fc_mx(
    eset,
    timepoint,
    timepoint_unit
  )

  # Get response
  response <- con$get_immune_response(
    assay = assay,
    participant_ids = rownames(FC),
    dichotomize = dichotomize,
    dichotomize_thresh = dichotomize_thresh
  )

  newdata <- data.frame(FC)
  rownames(newdata) <- rownames(FC)
  predicted_values <- stats::predict(fit, newdata = newdata, type = "response")

  predicted_values
}

#' Format Predictor Table
#'
#' Create nicely formatted table of predictors from relasso result.
#'
#' @param relasso fit object from \code{glm} or \code{lm}
#'
#' @export
format_predictor_table <- function(relasso) {
  sum_relasso <- summary(relasso)
  sum_relasso_coef <- sum_relasso$coefficients
  pred_cIdx <- grep("value|Pr", colnames(sum_relasso_coef))
  predictor_table <- sum_relasso_coef[, pred_cIdx][-1, ]
  colnames(predictor_table) <- c("statistic", "p-value")


  predictor_table <- data.table(cbind(
    data.table(selected_features = rownames(predictor_table)),
    predictor_table
  ))
  predictor_table[, .(
    gene_symbol = paste0(
      '<a href="http://immunet.princeton.edu/predictions/gene/?network=immune_global&gene=',
      selected_features,
      '" target="_blank">',
      selected_features,
      "</a>"
    ),
    statistic = as.numeric(statistic),
    `p-value` = as.numeric(`p-value`)
  )]
}
