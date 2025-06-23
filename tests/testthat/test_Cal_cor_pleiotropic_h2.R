test_that("Cal_cor_pleiotropic_h2 computes correct pleioh2g", {
  data(Results_full_rg_15D)
  data(h2_vector_15D)

  result <- Cal_cor_pleiotropic_h2(Results_full_rg_15D, h2_vector_15D)

  expect_type(result, "double")
  expect_length(result, 15)

  expect_false(any(is.na(result)))
})
