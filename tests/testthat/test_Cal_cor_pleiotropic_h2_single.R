test_that("Cal_cor_pleiotropic_h2_single computes correct pleioh2g", {
  data(Results_full_rg_15D)
  data(h2_vector_15D)

  plei_h2_idx <- 1
  h2g_T_single<-h2_vector_15D[plei_h2_idx]

  result <- Cal_cor_pleiotropic_h2_single(Results_full_rg_15D, h2g_T_single, plei_h2_idx)


  expect_type(result, "double")

  expect_false(is.na(result))
})
