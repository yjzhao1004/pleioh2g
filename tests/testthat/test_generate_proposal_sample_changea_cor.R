test_that("generate_proposal_sample_changea_cor computes correct noisy_inversed_element", {
  data(Results_full_rg)
  data(Results_full_rg_array)
  Results_full_rg<-Results_full_rg[1:15,1:15]
  Results_full_rg_array<-Results_full_rg_array[1:15,1:15,]
  plei_h2_idx <- 1

  ratio_a <- 0.75

  result <- generate_proposal_sample_changea_cor(Results_full_rg, Results_full_rg_array, plei_h2_idx, ratio_a)

  expect_type(result, "double")

  expect_false(is.na(result))
})
