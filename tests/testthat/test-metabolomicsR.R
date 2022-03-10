test_that("Load the dataset", {
  
  file_path <- system.file("extdata", "QMDiab_metabolomics_OrigScale.xlsx", package = "metabolomicsR", mustWork = TRUE)
  df_plasma <- load_excel(path = file_path,
                          data_sheet = 1,
                          feature_sheet = 4,
                          sample_sheet = 8,
                          sampleID = "QMDiab-ID",
                          featureID = "BIOCHEMICAL"
  )
})
