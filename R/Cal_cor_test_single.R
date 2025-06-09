#' Compute inversed elements for the target disease in bias correction procedure with referred disease index in the rg matrix
#'
#' This function inversed elements for the target disease in bias correction procedure.
#'
#' @param rg_mat genetic correlation matrix.
#' @param plei_h2_idx index of the target disease in the rg_mat.
#'
#' @return inverse element value for the target disease used for bias correction
#' @export
#'
#' @examples
#' data(Results_full_rg_15D)
#' plei_h2_idx<-1
#' Cal_cor_test_single(Results_full_rg_15D,plei_h2_idx)

Cal_cor_test_single <- function(rg_mat,plei_h2_idx){
  n <- dim(rg_mat)[1]
  g2dg1 <- rg_mat[plei_h2_idx, ][-plei_h2_idx]
  inverse_element_test <- t(g2dg1) %*% solve(rg_mat[-plei_h2_idx, -plei_h2_idx]) %*% g2dg1
  return(inverse_element_test)
}
