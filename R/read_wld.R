#' Read wld from either internal or external file - refer to ldscr R package (https://github.com/mglev1n/ldscr)
#'
#' @description
#' `read_wld()` Read wld from either internal or external file
#' @param wld (character) Path to directory containing weight files. Default is `NA`, which will utilize the built-in weight files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @return
#' A data frame (tibble) containing LD weight information read from the specified directory.
#' Each row corresponds to a SNP, and columns typically include:
#' \itemize{
#'   \item \code{CHR}: Chromosome number.
#'   \item \code{SNP}: SNP identifier (rsID).
#'   \item \code{BP}: Base pair position.
#'   \item \code{wLD}: Weight for LD regression.
#' }
#' @export
read_wld <- function(wld) {

    w <- fs::dir_ls(wld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())

  return(w)
}
