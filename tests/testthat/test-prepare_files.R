.prep_geo_files <- HIPCMatrix:::.prep_geo_files

test_that("immport files are loaded correctly", {
  # SDY1630: RNA-seq
  study <- "SDY1630"
  gef <- readRDS("test_data/sdy1630/SDY1630_gef.rds")
  platform <- "NA"
  input_files <- "GSE133383_filtered_NKcounts_table.748979.csv"
  analysis_dir <- "test_data/sdy1630"
  files <- HIPCMatrix:::.prep_immport_files(
    study,
    gef,
    platform,
    file.path(analysis_dir, input_files),
    analysis_dir,
    verbose
  )
  expect_equal(files, "test_data/sdy1630/supp_files/Spleen_ARM5244/SDY1630_raw_expression.txt")

  # SDY112: Affy
  study <- "SDY112"
  gef <- readRDS("test_data/sdy112/SDY112_gef.rds")
  analysis_dir <- "test_data/sdy112"
  platform <- "Affymetrix"
  input_files <- HIPCMatrix:::.select_input_files(analysis_dir)

  # Should do nothing
  files <- HIPCMatrix:::.prep_immport_files(
    study,
    gef,
    platform,
    input_files,
    analysis_dir,
    verbose
  )
  expect_true(all(input_files == files))

  # should error when one of the files is not CEL file
  input_files <- list.files(analysis_dir)
  expect_error(files <- HIPCMatrix:::.prep_immport_files(
    study,
    gef,
    platform,
    input_files,
    analysis_dir,
    verbose
  ),
  "Affymetrix should only have CEL input files")

  # SDY645: Illumina
  study <- "SDY387"
  input_files <- "SDY387_EXP13732_microarray.703288.tsv"
  gef <- readRDS("test_data/sdy387/SDY387_gef.rds")
  analysis_dir <- "test_data/sdy387"
  platform <- "Illumina"
  input_files <- unique(gef$file_info_name)
  files <- HIPCMatrix:::.prep_immport_files(
    study,
    gef,
    platform,
    file.path(analysis_dir, input_files),
    analysis_dir,
    verbose = FALSE
  )
  expect_equal(files,
               "test_data/sdy387/supp_files/Whole blood_ARM2325/SDY387_raw_expression.txt")
})

test_that(".prep_geo_files GSM SOFT", {
  meta_data <- list(file_location = "gsm_soft",
                    platform = "Illumina",
                    gsm_table_var_name = "AVERAGE_SIGNAL")
  gef <- con_all$getDataset("gene_expression_files",
                        colFilter = Rlabkey:::makeFilter(c("biosample_accession", "IN", "BS974298")),
                        original_view = TRUE)
  analysis_dir <- "test_data/sdy1289"
  expect_message(input_files <- .prep_geo_files("SDY1289",
                                                             gef,
                                                             meta_data,
                                                             input_files = NA,
                                                             analysis_dir,
                                                             verbose = TRUE,
                                                             reload = FALSE),
                 "Downloading GSM soft files to test_data/sdy1289/supp_files/Whole blood_ARM4455")
  expect_equal(input_files, "test_data/sdy1289/supp_files/Whole blood_ARM4455/SDY1289_raw_expression.txt")
})
test_that(".prep_geo_files GSM supp illumina", {
  meta_data <- get_meta_data("SDY180")
  gef <- con_all$getDataset("gene_expression_files",
                            colFilter = Rlabkey:::makeFilter(c("biosample_accession", "IN", "BS662409")),
                            original_view = TRUE)[!is.na(geo_accession)]
  analysis_dir <- "test_data/sdy180"
  expect_message(input_files <- .prep_geo_files("SDY180",
                                              gef,
                                              meta_data,
                                              input_files = NA,
                                              analysis_dir,
                                              verbose = TRUE,
                                              reload = FALSE),
                 "Downloading GSM supp files to test_data/sdy180/supp_files/Whole blood_ARM773")
  expect_equal(input_files, "test_data/sdy180/supp_files/Whole blood_ARM773/SDY180_raw_expression.txt")
})
test_that(".prep_geo_files GSM supp affy", {
  meta_data <- get_meta_data("SDY1328")
  gef <- con_all$getDataset("gene_expression_files",
                            colFilter = Rlabkey:::makeFilter(c("biosample_accession", "IN", "BS978363;BS1005596")),
                            original_view = TRUE)[!is.na(geo_accession)]
  analysis_dir <- "test_data/sdy1328"
  expect_message(input_files <- .prep_geo_files("SDY1328",
                                                gef,
                                                meta_data,
                                                input_files = NA,
                                                analysis_dir,
                                                verbose = TRUE,
                                                reload = FALSE),
                 "Downloading GSM supp files to test_data/sdy1328/supp_files/Whole blood_ARM4537")
  expect_length(input_files, 2)
  expect_true(all(grepl("CEL", input_files)))
})
test_that(".prep_geo_files GSM supp rnaseq", {
  meta_data <- list(file_location = "gsm_supp_files",
                    platform = "NA")
  gef <- con_all$getDataset("gene_expression_files",
                            colFilter = Rlabkey:::makeFilter(c("biosample_accession", "IN", "BS1004262;BS1004309")),
                            original_view = TRUE)[!is.na(geo_accession)]
  analysis_dir <- "test_data/sdy1412"
  expect_message(input_files <- .prep_geo_files("SDY1412",
                                                gef,
                                                meta_data,
                                                input_files = NA,
                                                analysis_dir,
                                                verbose = TRUE,
                                                reload = FALSE),
                 "GSM3494630")
  expect_true(grepl("SDY1412_raw_expression.txt", input_files))
  ge <- fread(input_files)
  expect_equal(ncol(ge), 3)
})


test_that(".prep_geo_files GSE supp illumina", {
  meta_data <- get_meta_data("SDY640")
  gef <- con_all$getDataset("gene_expression_files",
                            colFilter = Rlabkey:::makeFilter(c("biosample_accession", "IN", "BS799828")),
                            original_view = TRUE)[!is.na(geo_accession)]
  analysis_dir <- "test_data/sdy640"
  expect_message(input_files <- .prep_geo_files("SDY640",
                                 gef,
                                 meta_data,
                                 input_files = NA,
                                 analysis_dir,
                                 verbose = TRUE,
                                 reload = FALSE),
                 "Using saved supp files for GSE101710")
  expect_true(grepl("SDY640_raw_expression.txt", input_files))
})
test_that(".prep_geo_files GSE supp rnaseq", {
  meta_data <- get_meta_data("SDY787")
  gef <- readRDS("test_data/sdy787/SDY787_gef.rds")[1, ]
  analysis_dir <- "test_data/sdy787"
  expect_message(input_files <- HIPCMatrix:::.prep_geo_files("SDY787",
                                              gef,
                                              meta_data,
                                              input_files = NA,
                                              analysis_dir,
                                              verbose = TRUE,
                                              reload = FALSE),
                 "Using locally cached version of GSM3122900")
  expect_equal(input_files, "test_data/sdy787/supp_files/T cell_ARM3144/SDY787_raw_expression.txt")
})


get_platform <- function(study) {
  unique(Rlabkey::labkey.selectRows(baseUrl = labkey.url.base,
                                    folderPath = paste0("Studies/", study),
                                    schemaName = "assay.ExpressionMatrix.matrix",
                                    queryName = "Runs",
                                    colNameOpt = "fieldname",
                                    colSelect = "featureSet/Vendor")$`featureSet/Vendor`)
}

test_that("select_input_files selects the correct files", {
  supp_files_dir <- "test_data/sdy787/supp_files/T cell_ARM3144"
  expect_equal(HIPCMatrix:::.select_input_files(supp_files_dir),
               file.path(supp_files_dir, "GSE113891_PT_all_Count.tsv"))
})
