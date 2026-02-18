#' Internal function to perform LDSC heritability/covariance analysis - refer to ldscr R package (https://github.com/mglev1n/ldscr)

#' @description
#' `perform_analysis()` Internal function to perform LDSC heritability/covariance analysis

#' @param n.blocks Number of blocks
#' @param n.snps Number of SNPs
#' @param weighted.LD wld score
#' @param weighted.chi chi-square
#' @param N.bar Average N after merging
#' @param m Number of SNPs from LD data
#' @return
#' A list containing the results of the LDSC heritability/covariance analysis with the following elements:
#' \itemize{
#'   \item \code{reg.tot}: Estimated total heritability or covariance (regression coefficient scaled by \code{m}).
#'   \item \code{tot.se}: Standard error of the total heritability/covariance estimate, computed using a block jackknife.
#'   \item \code{intercept}: LDSC regression intercept.
#'   \item \code{intercept.se}: Standard error of the intercept, estimated via block jackknife.
#'   \item \code{pseudo.values}: Vector of pseudo-values from the block jackknife procedure, one per block.
#'   \item \code{N.bar}: Average sample size across SNPs after merging.
#' }
#' @export

perform_analysis <- function(n.blocks, n.snps, weighted.LD, weighted.chi, N.bar, m) {
  n.annot <- 1

  select.from <- floor(seq(from = 1, to = n.snps, length.out = (n.blocks + 1)))
  select.to <- c(select.from[2:n.blocks] - 1, n.snps)

  xty.block.values <- matrix(data = NA, nrow = n.blocks, ncol = (n.annot + 1))
  xtx.block.values <- matrix(data = NA, nrow = ((n.annot + 1) * n.blocks), ncol = (n.annot + 1))
  colnames(xty.block.values) <- colnames(xtx.block.values) <- colnames(weighted.LD)
  replace.from <- seq(from = 1, to = nrow(xtx.block.values), by = (n.annot + 1))
  replace.to <- seq(from = (n.annot + 1), to = nrow(xtx.block.values), by = (n.annot + 1))
  for (i in 1:n.blocks) {
    xty.block.values[i, ] <- t(t(weighted.LD[select.from[i]:select.to[i], ]) %*% weighted.chi[select.from[i]:select.to[i], ])
    xtx.block.values[replace.from[i]:replace.to[i], ] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i], ]) %*% weighted.LD[select.from[i]:select.to[i], ])
  }
  xty <- as.matrix(colSums(xty.block.values))
  xtx <- matrix(data = NA, nrow = (n.annot + 1), ncol = (n.annot + 1))
  colnames(xtx) <- colnames(weighted.LD)
  for (i in 1:nrow(xtx)) {
    xtx[i, ] <- t(colSums(xtx.block.values[seq(from = i, to = nrow(xtx.block.values), by = ncol(weighted.LD)), ]))
  }

  reg <- solve(xtx, xty)
  intercept <- reg[2]
  coefs <- reg[1] / N.bar
  reg.tot <- coefs * m

  delete.from <- seq(from = 1, to = nrow(xtx.block.values), by = ncol(xtx.block.values))
  delete.to <- seq(from = ncol(xtx.block.values), to = nrow(xtx.block.values), by = ncol(xtx.block.values))
  delete.values <- matrix(data = NA, nrow = n.blocks, ncol = (n.annot + 1))
  colnames(delete.values) <- colnames(weighted.LD)
  for (i in 1:n.blocks) {
    xty.delete <- xty - xty.block.values[i, ]
    xtx.delete <- xtx - xtx.block.values[delete.from[i]:delete.to[i], ]
    delete.values[i, ] <- solve(xtx.delete, xty.delete)
  }

  tot.delete.values <- delete.values[, 1:n.annot]
  pseudo.values <- matrix(data = NA, nrow = n.blocks, ncol = length(reg))
  colnames(pseudo.values) <- colnames(weighted.LD)
  for (i in 1:n.blocks) {
    pseudo.values[i, ] <- (n.blocks * reg) - ((n.blocks - 1) * delete.values[i, ])
  }

  jackknife.cov <- cov(pseudo.values) / n.blocks
  jackknife.se <- sqrt(diag(jackknife.cov))
  intercept.se <- jackknife.se[length(jackknife.se)]
  coef.cov <- jackknife.cov[1:n.annot, 1:n.annot] / (N.bar^2)
  cat.cov <- coef.cov * (m %*% t(m))
  tot.cov <- sum(cat.cov)
  tot.se <- sqrt(tot.cov)

  return(
    list(
      reg.tot = reg.tot,
      tot.se = tot.se,
      intercept = intercept,
      intercept.se = intercept.se,
      pseudo.values = pseudo.values[, 1],
      N.bar = N.bar
    )
  )
}
