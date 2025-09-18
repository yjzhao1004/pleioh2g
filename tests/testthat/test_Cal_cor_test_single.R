test_that("Cal_cor_test_single computes correct inverse element", {
  data(Results_full_rg_15D)

  plei_h2_idx <- 1

  result <- Cal_cor_test_single(Results_full_rg_15D, plei_h2_idx)

  expect_type(result, "double")

  expect_false(is.na(result))
})
