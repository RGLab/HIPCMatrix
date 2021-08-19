#' HIPCMatrix Connection
#'
#' @description ImmuneSpace connection object with additional methods
#' for gene expression analysis
#'
#' @importFrom R6 R6Class
#' @export
HMX <- R6::R6Class(
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
    }
  )
)
