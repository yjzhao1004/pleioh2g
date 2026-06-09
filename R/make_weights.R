#' Internal Function to make weights - refer to ldscr R package (https://github.com/mglev1n/ldscr)

#' @description
#' `make_weights()` Internal Function to make weights

#' @param chi1 chi-square
#' @param L2 ld score
#' @param wLD wld score
#' @param N sample size
#' @param M.tot Number of SNPs
#' @return
#' A numeric vector of initial LDSC weights for each SNP
#'
#' @export
#'
make_weights <- function(chi1, L2, wLD, N, M.tot) {
  tot.agg <- (M.tot * (mean(chi1) - 1)) / mean(L2 * N)
  tot.agg <- max(tot.agg, 0)
  tot.agg <- min(tot.agg, 1)
  ld <- pmax(L2, 1)
  w.ld <- pmax(wLD, 1)
  c <- tot.agg * N / M.tot
  het.w <- 1 / (2 * (1 + (c * ld))^2)
  oc.w <- 1 / w.ld
  w <- het.w * oc.w
  initial.w <- sqrt(w)

  return(initial.w)
}
