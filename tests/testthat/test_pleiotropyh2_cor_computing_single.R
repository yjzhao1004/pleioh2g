test_that("pleiotropyh2_cor_computing_single works correctly", {

  library(testthat)
  library(data.table)

  data(Results_full_rg_15D)
  data(Results_full_rg_array_15D)
  data(h2_vector_15D)
  data(h2_vector_mat_15D)
  num_diseases <- ncol(h2_vector_15D)
  phenotype_path <- system.file("extdata", "phenotype_name_package.txt", package = "pleioh2g")
  save_path <- system.file("extdata", "save_results_test/D",package = "pleioh2g")

  pleiotropyh2_cor_computing_single_prune(1,phenotype_path,save_path,h2_vector_15D,h2_vector_mat_15D,Results_full_rg_15D,Results_full_rg_array_15D, sample_rep = 10, tolerance = 1e-6,seed = 123)


  result_rds <- file.path(save_path, "401.1_perpleioh2g.rds")
  expect_true(file.exists(result_rds))

  result_file <- file.path(save_path, "401.1_h2pleio_dresults.csv")
  expect_true(file.exists(result_file))

  results <- read.csv(result_file)
  expect_true("target_disease" %in% colnames(results))
  expect_true("target_disease_h2_est" %in% colnames(results))
  expect_true("target_disease_h2_se" %in% colnames(results))
  expect_true("selected_pleio_pheno" %in% colnames(results))
  expect_true("h2pleiotropy" %in% colnames(results))
  expect_true("h2pleiotropy_se" %in% colnames(results))
  expect_true("percentage_jackknife" %in% colnames(results))
  expect_true("percentage_jackknife_se" %in% colnames(results))
  expect_true("h2pleio_corr" %in% colnames(results))
  expect_true("percentage_corr" %in% colnames(results))
  expect_true("corrected_weight" %in% colnames(results))


})
