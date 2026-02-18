test_that("Cal_cor_pleiotropic_h2_corrected_single computes correct pleioh2g after correction", {
  data(Results_full_rg)
  data(h2_vector)
  plei_h2_idx <- 1

  h2g_T_single <- h2_vector[plei_h2_idx]

  corrected_weight_updated <- 0.78

  result <- Cal_cor_pleiotropic_h2_corrected_single(Results_full_rg, h2g_T_single, corrected_weight_updated, plei_h2_idx)

  expect_type(result, "double")

  expect_false(is.na(result))

})
