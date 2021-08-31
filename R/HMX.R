#' HIPCMatrix Connection
#'
#' @description ImmuneSpace connection object with additional methods
#' for gene expression analysis
#'
#' @importFrom R6 R6Class
#' @export
HMX <- R6Class(
  classname = "HIPCMatrixConn",
  inherit = ImmuneSpaceR:::ISCon,
  public = list(

    # Methods


    #' @description Run Differential Expression analysis on matrices in a study
    #' @aliases runGEAnalysis
    #'
    #' @param rerun Force re-run if results already in cache?
    #'
    runGEAnalysis = function(rerun = FALSE) {
      runGEAnalysis(self, rerun)
    },


    #' @description Upload differential gene expression analysis results to
    #' server. This updates the gene_expression_analysis and the
    #' gene_expression_analysis_results tables with updated results.
    #' @aliases uploadGEAnalysisResults
    uploadGEAnalysisResults = function() {
      uploadGEAnalysisResults(self)
    },

    #' @description Check gene expression matrices to determine if
    #' differential expression can be run.
    #'
    #' @return `TRUE` if differential expression can be run.
    #' `FALSE` if not.
    checkImpliedGEAR = function() {
      checkImpliedGEAR(self)
    },

    #' @description Update Microarray.FeatureAnnotation table
    #' @aliases updateFAS
    #'
    #' @param fas_names name of feature annotation set
    #' @details  This will take care of GeneExpressionExplorer Module, which uses
    #' the Microarray.FeatureAnnotation table to populate a dropdown for selection
    #' of genes of interest.  GEE will look for the FASid that was given to the original
    #' FAS when it was uploaded and this is challenging to change, therefore
    #' it was decided to move the original to a new FAS called "myOriginalFasName_orig"
    #' that gets a new FASid.  Con$getGEMatrix() then looks for this new FASid when
    #' populating the probe level gene symbols with the arg: `annotation = "default"`.
    updateFAS = function(fas_names) {
      updateFAS(self, fas_names)
    },

    #' @description Update Expression Matrices with new annotation
    #' @aliases updateEMs
    #'
    #' @details This will update the summary.tsv files for all matrices
    #' associated with the connection by re-calculating the summarized
    #' expression values with the most current prob to gene symbol mapping.
    updateEMs = function() {
      updateEMs(self)
    },

    #' @description Update gene annotation
    #' @aliases runUpdateAnno
    #'
    #' @details This function will update the feature annotation tables,
    #' summarized expression matrices, and re-run differential expression
    #' where relevant, using the most current gene alias to symbol mapping
    #' found in the HIPCMatrix package data.
    #'
    #' To update annotation across ImmuneSpace, first update and re-install
    #' UpdateAnno. Then:
    #'
    #' \code{con <- HMX$new("")}
    #'
    #' \code{con$runUpdateAnno()}
    #'
    #' To update feature annotation and run differential expression on
    #' a new matrix after loading all matrices for a study:
    #'
    #' \code{con <- HMX$new("SDY_ID")}
    #'
    #' \code{con$runUpdateAnno()}
    runUpdateAnno = function() {
      if (grepl("^IS", self$study)) {
        stop("Do not update ImmuneSignatures Annotation!")
      }

      self$updateFAS()

      if (self$study == "Studies") {
        ge_studies <- unique(self$listGEMatrices()$folder)

        lapply(ge_studies, function(study) {
          con <- HMX$new(
            study = study,
            verbose = TRUE
          )
          con$updateEMs()
          con$uploadGEAnalysisResults()
        })
      } else {
        self$updateEMs()
        self$uploadGEAnalysisResults()
      }
      invisible(self)
    },


    #' @description Get immune response.
    #' @aliases get_immune_response
    #'
    #' @return Max log fold change of \code{assay} from baseline for each
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
    #' @param reload Force rerun if result is already found in cache?
    get_immune_response = function(assay = "hai",
                                   participant_ids = NULL,
                                   dichotomize = FALSE,
                                   dichotomize_thresh = 4,
                                   reload = FALSE) {
      get_immune_response(
        self,
        assay = assay,
        participant_ids = participant_ids,
        dichotomize = dichotomize,
        dichotomize_thresh = dichotomize_thresh,
        reload = reload
      )
    },


    #' @description Get genes differentially expressed over time
    #' @aliases get_de_genes
    #'
    #' @param timepoint Timepoint for finding differentially expressed genes
    #' @param fc_thresh fold-change threshold for determining whether genes
    #' are differentially expressed
    #' @param timepoint_unit "Days" or "Hours"
    #' @param cohorts character vector of cohorts to include
    #' @param reload Force rerun if result is already found in cache?
    get_de_genes = function(timepoint,
                            fc_thresh = 0.58,
                            timepoint_unit = "Days",
                            cohorts = NULL,
                            reload = FALSE) {
      get_de_genes(
        self,
        timepoint = timepoint,
        fc_thresh = fc_thresh,
        timepoint_unit = timepoint_unit,
        cohorts = cohorts,
        reload = reload
      )
    },

    #' @description Get immune response predictors
    #' @aliases train_immune_response_predictors
    #'
    #' @param cohorts character vector of cohorts to include
    #' @param timepoint Timepoint for finding differentially expressed genes
    #' @param assay "hai", "neut_ab_response", or "elisa"
    #' @param timepoint_unit "Days" or "Hours"
    #' @param use_only_de_genes Filter to differentially expressed genes
    #' when finding predictos?
    #' @param fc_thresh fold-change threshold for determining whether genes
    #' are differentially expressed. Ignored if \code{use_only_de_genes} is
    #' \code{FALSE}.
    #' @param dichotomize Dichotomize result? If FALSE, max log fold change is
    #' returned. If TRUE, returns TRUE if max log fold change is greater than
    #' \code{dichotomize_value}
    #' @param dichotomize_thresh Value to use for dichotomizing result if
    #' \code{dichotomize} is TRUE.
    #' @param return_type Option to return either a vector of the significant
    #' features, or the model fit. Options are 'fit' or 'features'. Default
    #' is to return the fitted model.
    #' @param reload Force rerun if result is already found in cache?
    train_immune_response_predictors = function(cohorts,
                                                timepoint,
                                                assay = "hai",
                                                timepoint_unit = "Days",
                                                use_only_de_genes = timepoint > 0,
                                                fc_thresh = 0.58,
                                                dichotomize = FALSE,
                                                dichotomize_thresh = 4,
                                                return_type = "fit",
                                                reload = FALSE) {
      train_immune_response_predictors(
        self,
        cohorts = cohorts,
        timepoint = timepoint,
        assay = assay,
        timepoint_unit = timepoint_unit,
        use_only_de_genes = use_only_de_genes,
        fc_thresh = fc_thresh,
        dichotomize = dichotomize,
        dichotomize_thresh = dichotomize_thresh,
        return_type = return_type,
        reload = reload
      )
    },

    #' @description Test immune response predictor model on testing cohort data.
    #' @aliases predict_response
    #'
    #' @param cohorts character vector of cohorts to include
    #' @param timepoint Timepoint for finding differentially expressed genes
    #' @param fit Model fit to test.
    #' @param timepoint_unit "Days" or "Hours"
    predict_response = function(cohorts,
                                timepoint,
                                fit,
                                timepoint_unit = "Days") {
      predict_response(self,
        cohorts,
        timepoint,
        fit,
        timepoint_unit = "Days"
      )
    },

    #' @description Test immune response predictors
    #'
    #' @details get a table of observed vs predicted
    #' values given a fitted model.
    #'
    #' @param cohorts cohorts to test
    #' @param timepoint Timepoint to test
    #' @param fit Model fit to test
    #' @param assay "hai", "neut_ab_response", or "elisa"
    #' @param timepoin_unit "Days" or "Hours"
    #' @param dichotomize Dichotomize result? If FALSE, max log fold change is
    #' returned. If TRUE, returns TRUE if max log fold change is greater than
    #' \code{dichotomize_value}
    #' @param dichotomize_thresh Value to use for dichotomizing result if
    #' \code{dichotomize} is TRUE.
    test_immune_response_predictors = function(cohorts,
                                               timepoint,
                                               fit,
                                               assay = "hai",
                                               timepoint_unit = "Days",
                                               dichotomize = FALSE,
                                               dichotomize_thresh = 4) {
      test_immune_response_predictors(
        self,
        cohorts = cohorts,
        timepoint = timepoint,
        fit = fit,
        assay = assay,
        timepoint_unit = timepoint_unit,
        dichotomize = dichotomize,
        dichotomize_thresh = dichotomize_thresh
      )
    }
  )
)
