#' Read ld from either internal or external file - refer to ldscr R package (https://github.com/mglev1n/ldscr)
#'
#' @description
#' `read_ld()` Read ld from either internal or external file.
#' @param ld (character) Path to directory containing ld score files, ending in `*.l2.ldscore.gz`. Default is `NA`, which will utilize the built-in ld score files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @return
#' A data frame (tibble) containing LD score information read from the specified directory.
#' Each row corresponds to a SNP, and columns typically include:
#' \itemize{
#'   \item \code{CHR}: Chromosome number.
#'   \item \code{SNP}: SNP identifier (rsID).
#'   \item \code{BP}: Base pair position.
#'   \item \code{L2}: LD score value.
#'   \item \code{M}: Number of SNPs used in the LD score computation.
#' }
#'
#' @export
read_ld <- function(ld) {

    x <- fs::dir_ls(ld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())

  return(x)
}
