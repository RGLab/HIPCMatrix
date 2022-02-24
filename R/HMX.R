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

    # Properties


    #' @field immune_response_predictors List of ImmuneResponsePredictor
    #' objects associated with this connection. Created by \code{HMS$run_irp()}
    immune_response_predictors = list(),

    # Methods


    #' @description Run Differential Expression analysis on matrices in a study
    #' @aliases run_de_analysis
    #'
    #' @param rerun Force re-run if results already in cache?
    #'
    run_de_analysis = function(rerun = FALSE) {
      run_de_analysis(self, rerun)
    },


    #' @description Upload differential gene expression analysis results to
    #' server. This updates the gene_expression_analysis and the
    #' gene_expression_analysis_results tables with updated results.
    #' @aliases uploadGEAnalysisResults
    upload_de_analysis_results = function() {
      upload_de_analysis_results(self)
    },

    #' @description Find gene expression matrices and timepoints which are compatible
    #' with differential expression analysis
    #'
    #' @return data.table with one row per timepoint per matrix which has
    #' sufficient data to run differential expression analysis
    get_de_compatible_runs = function() {
      get_de_compatible_runs(self)
    },

    #' @description Update Microarray.FeatureAnnotation table
    #' @aliases updateFAS
    #'
    #' @param fas_names name of feature annotation set. If NULL, any FAS
    #' associated with matrices associated with the connection will be updated.
    #' For global connection, all FAS will be updated.
    #' @details  This will take care of GeneExpressionExplorer Module, which uses
    #' the Microarray.FeatureAnnotation table to populate a dropdown for selection
    #' of genes of interest.  GEE will look for the FASid that was given to the original
    #' FAS when it was uploaded and this is challenging to change, therefore
    #' it was decided to move the original to a new FAS called "myOriginalFasName_orig"
    #' that gets a new FASid.  Con$getGEMatrix() then looks for this new FASid when
    #' populating the probe level gene symbols with the arg: `annotation = "default"`.
    updateFAS = function(fas_names = NULL) {
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


    #' @description Get immune response.
    #' @aliases get_immune_response
    #'
    #' @return Max log fold change of \code{assay} from baseline for each
    #' particiapnt.
    #'
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
    #' @return Invisibly returns \code{irp_index} which can be used to
    #' access the ImmuneResponsePredictor object with \code{HMX$get_irp()}
    #'
    #' @param cohorts character vector of cohorts to include
    #' @param timepoint Timepoint for finding differentially expressed genes
    #' @param assay "hai", "neut_ab_response", or "elisa"
    #' @param timepoint_unit "Days" or "Hours"
    #' @param use_only_de_genes Filter to differentially expressed genes
    #' when finding predictors?
    #' @param fc_thresh fold-change threshold for determining whether genes
    #' are differentially expressed. Ignored if \code{use_only_de_genes} is
    #' \code{FALSE}.
    #' @param dichotomize Dichotomize result? If FALSE, max log fold change is
    #' returned. If TRUE, returns TRUE if max log fold change is greater than
    #' \code{dichotomize_value}
    #' @param dichotomize_thresh Value to use for dichotomizing result if
    #' \code{dichotomize} is TRUE.
    #' @param reload Force rerun if result is already found in cache?
    train_immune_response_predictors = function(cohorts,
                                                timepoint,
                                                assay = "hai",
                                                timepoint_unit = "Days",
                                                use_only_de_genes = timepoint > 0,
                                                fc_thresh = 0.58,
                                                dichotomize = FALSE,
                                                dichotomize_thresh = 4,
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
        reload = reload
      )
    },

    #' @description Get ImmuneResponsePredictor
    #'
    #' @param irp_index index of immune_response_predictor object in HMX
    get_irp = function(irp_index = NULL) {
      if (is.null(irp_index)) irp_index <- length(self$immune_response_predictors)
      self$immune_response_predictors[[irp_index]]
    },

    #' @description Test immune response predictor model on testing cohort data.
    #' @aliases predict_response
    #'
    #' @param cohort cohort to use (string)
    #' @param irp_index index of immune_response_predictor object in HMX
    predict_response = function(cohort,
                                irp_index = NULL) {
      predict_response(
        self,
        cohort,
        irp_index
      )
    },

    #' @description get a table of observed vs predicted
    #' values given a fitted model.
    #'
    #' @param cohorts cohorts to test
    #' @param irp_index irp_index
    #' @param reload force re-run if already cached?
    test_immune_response_predictors = function(cohorts,
                                               irp_index = NULL,
                                               reload = FALSE) {
      test_immune_response_predictors(
        self,
        cohorts = cohorts,
        irp_index = irp_index,
        reload = reload
      )
    },

    #' @description Find predictors of immune response from gene expression.
    #' @aliases run_irp
    #'
    #' @param cohorts_train Training cohorts
    #' @param cohorts_test Testing cohorts (optional)
    #' @param timepoint Timepoint for finding differentially expressed genes
    #' @param assay "hai", "neut_ab_response", or "elisa"
    #' @param timepoint_unit "Days" or "Hours"
    #' @param use_only_de_genes Filter to differentially expressed genes
    #' when finding predictors?
    #' @param fc_thresh fold-change threshold for determining whether genes
    #' are differentially expressed. Ignored if \code{use_only_de_genes} is
    #' \code{FALSE}.
    #' @param dichotomize Dichotomize result? If FALSE, max log fold change is
    #' returned. If TRUE, returns TRUE if max log fold change is greater than
    #' \code{dichotomize_value}
    #' @param dichotomize_thresh Value to use for dichotomizing result if
    #' \code{dichotomize} is TRUE.
    #' @param reload Force rerun if result is already found in cache?
    run_irp = function(cohorts_train,
                       cohorts_test = NULL,
                       timepoint,
                       use_only_de_genes = timepoint > 0,
                       assay = "hai",
                       timepoint_unit = "Days",
                       fc_thresh = 0.58,
                       dichotomize = FALSE,
                       dichotomize_thresh = 4,
                       reload = FALSE) {
      run_irp(
        self,
        cohorts_train = cohorts_train,
        cohorts_test = cohorts_test,
        timepoint = timepoint,
        use_only_de_genes = use_only_de_genes,
        assay = assay,
        timepoint_unit = timepoint_unit,
        fc_thresh = fc_thresh,
        dichotomize = dichotomize,
        dichotomize_thresh = dichotomize_thresh,
        reload = reload
      )
    },

    #' @description Run a gene set enrichment analysis on a gene expression matrix,
    #' comparing all timepoints to baseline. The CAMERA method from the limma
    #' package is used, as described in
    #' \url{https://academic.oup.com/nar/article/40/17/e133/2411151}{Wu and Smyth (2012)}
    #'
    #' @aliases run_gsea
    #'
    #'
    #' @param matrix_name The name of the gene expression matrix
    #'  to download.
    #' @param cohort_type  The name of a cohortType that has an
    #' associated gene expression matrix. Note that if this argument is not
    #' NULL, then \code{matrixName} is ignored. CohortType is a concatenation of
    #' "cohort" and "cell type" that allows the user to specify a matrix for the
    #' cell type subset of a cohort.
    #' @param set_name Name of predefined set of gene signatures. Choose from:
    #' \code{chaussabel}, \code{blood_transcription}, \code{msigdb}
    #' @param gene_sets  A list of vectors of gene names, each entry corresponding
    #' to a gene set. If specified, this will be used in place of the "set_name"
    #' argument to test gene sets.
    run_gsea = function(matrix_name = NULL,
                        cohort_type = NULL,
                        set_name = "msigdb",
                        gene_sets = NULL) {
      run_gsea(
        con = self,
        matrix_name = matrix_name,
        cohort_type = cohort_type,
        set_name = set_name,
        gene_sets = gene_sets
      )
    }
  )
)
