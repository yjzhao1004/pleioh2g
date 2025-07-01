#' Compute a vector of pleioh2g for all diseases before correction
#' This function computes pleioh2g for all diseases before correction in one go.
#' @param rg_mat genetic correlation matrix.
#' @param h2g_T heritability vector for all diseases.
#'
#' @return pleioh2g vector
#' @export
#'
#' @examples
#' data(Results_full_rg_15D)
#' data(h2_vector_15D)
#' Cal_cor_pleiotropic_h2(Results_full_rg_15D,h2_vector_15D)

Cal_cor_pleiotropic_h2 <- function(rg_mat,h2g_T){
  n <- dim(rg_mat)[1]
  h2g <- rep(NA, n)
  for(plei_h2_idx in 1:n){ # pick one disease to compute pleiotropic h2g
    g2dg1 <- rg_mat[plei_h2_idx, ][-plei_h2_idx]
    cov_g2d <- rg_mat[-plei_h2_idx, -plei_h2_idx]
    h2g[plei_h2_idx] <- h2g_T[plei_h2_idx]*(t(g2dg1) %*% solve(rg_mat[-plei_h2_idx, -plei_h2_idx]) %*% g2dg1)
  }
  return(h2g)
}
