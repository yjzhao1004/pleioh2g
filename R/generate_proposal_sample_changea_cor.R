#' Generate samples based on sampling covariance matrix and rg matrix for target disease
#'
#' This function is used to generate samples based on sampling covariance matrix and rg matrix for target disease
#'
#' @param Results_full_rg genetic correlation matrix.
#' @param Results_full_rg_array genetic correlation jackknife-block array.
#' @param plei_h2_idx index of the target disease in the rg_mat.
#' @param ratio_a corrected ratio.
#'
#' @import mvtnorm
#' @importFrom stats cov
#' @importFrom stats sd
#' @importFrom utils write.csv

#' @return noisy_inversed_element for bias correction
#' @export
#'
#' @examples
#' data(Results_full_rg_15D)
#' data(Results_full_rg_array_15D)
#' plei_h2_idx<-1
#' ratio_a <- 0.75
#' generate_proposal_sample_changea_cor(Results_full_rg_15D, Results_full_rg_array_15D, plei_h2_idx, ratio_a)

generate_proposal_sample_changea_cor <- function(Results_full_rg, Results_full_rg_array, plei_h2_idx, ratio_a) {
  n = nrow(Results_full_rg)
  n_block<-dim(Results_full_rg_array)[3]
  cov_matrix <- matrix(0, nrow = n*(n-1)/2, ncol = n_block)
  for(i in c(1:n_block)){
    current_mat<-Results_full_rg_array[,,i]
    lower_tri_indices <- which(lower.tri(current_mat, diag = F), arr.ind = TRUE)
    lower_tri_elements <- current_mat[lower_tri_indices]
    cov_matrix[,i]<-lower_tri_elements
  }
  covcov<-cov(t(cov_matrix))
  covcov_scale<-covcov*n_block
  iter=1
  repeat {

    lower_noise_matrix <- matrix(0, n, n)
    lower_noise_matrix[lower.tri(lower_noise_matrix, diag = F)] <- mvtnorm::rmvnorm(1, sigma = covcov_scale)
    upper_triangle_noise_matrix <- t(lower_noise_matrix)
    noisy_matrix<-lower_noise_matrix + upper_triangle_noise_matrix
    sample_rg_noise <- Results_full_rg + noisy_matrix
    if (iter > 50) {

      return(NA)
      break
    }

    if (all(eigen(sample_rg_noise)$values >0) && all(diag(sample_rg_noise)>0)) {
      sample_a_noise <- ratio_a*Results_full_rg[plei_h2_idx, ][-plei_h2_idx] + noisy_matrix[plei_h2_idx, ][-plei_h2_idx]
      sample_rgsolve_noise <- Results_full_rg[-plei_h2_idx, -plei_h2_idx] + noisy_matrix[-plei_h2_idx, -plei_h2_idx]
      noisy_inversed_element <- t(sample_a_noise) %*% solve(sample_rgsolve_noise) %*% sample_a_noise
      return(noisy_inversed_element)
    }
  }
}

