
#' Compute single pleioh2g for target disease after correction with referred disease index in the rg matrix and corrected ratio
#'
#' This function computes pleioh2g for the target disease after correction.
#'
#' @param rg_mat genetic correlation matrix.
#' @param h2g_T_single heritability for target diseases.
#' @param corrected_weight_updated the ratio for correction
#' @param plei_h2_idx index of the target disease in the rg_mat.
#'
#' @return pleioh2g value for the target disease after correction
#' @export
#'
#' @examples
#' data(Results_full_rg_15D)
#' data(h2_vector_15D)
#' plei_h2_idx<-1
#' h2g_T_single <- h2_vector_15D[plei_h2_idx]
#' corrected_weight_updated <- 0.78
#' Cal_cor_pleiotropic_h2_corrected_single(Results_full_rg_15D,h2g_T_single,
#' corrected_weight_updated,plei_h2_idx)

Cal_cor_pleiotropic_h2_corrected_single <- function(rg_mat,h2g_T_single,corrected_weight_updated,plei_h2_idx){
  g2dg1 <- rg_mat[plei_h2_idx, ][-plei_h2_idx]
  cov_g2d <- rg_mat[-plei_h2_idx, -plei_h2_idx]
  h2g <- h2g_T_single*((corrected_weight_updated*t(g2dg1)) %*% solve(rg_mat[-plei_h2_idx, -plei_h2_idx]) %*% (corrected_weight_updated*g2dg1))
  return(h2g)
}
