test_that("Cal_cor_pleiotropic_h2 computes correct pleioh2g", {
  data(Results_full_rg)
  data(h2_vector)

  result <- Cal_cor_pleiotropic_h2(Results_full_rg, h2_vector)

  expect_type(result, "double")
  expect_length(result, 62)

  expect_false(any(is.na(result)))
})
