
#' Compute single pleioh2g for target disease before correction with referred disease index in the rg matrix
#'
#' This function computes pleioh2g for the target disease before correction.
#'
#' @param rg_mat genetic correlation matrix.
#' @param h2g_T_single heritability for target diseases.
#' @param plei_h2_idx index of the target disease in the rg_mat.
#'
#' @return pleioh2g value for the target disease before correction
#' @export
#'
#' @examples
#' data(Results_full_rg_15D)
#' data(h2_vector_15D)
#' plei_h2_idx<-1
#' h2g_T_single<-h2_vector_15D[plei_h2_idx]
#' Cal_cor_pleiotropic_h2_single(Results_full_rg_15D,h2g_T_single,plei_h2_idx)

Cal_cor_pleiotropic_h2_single <- function(rg_mat,h2g_T_single,plei_h2_idx){
  n <- dim(rg_mat)[1]
  g2dg1 <- rg_mat[plei_h2_idx, ][-plei_h2_idx]
  cov_g2d <- rg_mat[-plei_h2_idx, -plei_h2_idx]
  h2g <- h2g_T_single*(t(g2dg1) %*% solve(rg_mat[-plei_h2_idx, -plei_h2_idx]) %*% g2dg1)
  return(h2g)
}
