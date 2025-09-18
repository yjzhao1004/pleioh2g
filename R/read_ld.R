#' Read ld from either internal or external file
#'
#' @description
#' `read_ld()` Read ld from either internal or external file.

#' @param ancestry (character) One of "AFR", "AMR", "CSA", "EAS", "EUR", or "MID", which will utilize the appropriate built-in `ld` and `wld` files from Pan-UK Biobank. If empty or `NULL`, the user must specify paths to `ld` and `wld` files.
#' @param ld (character) Path to directory containing ld score files, ending in `*.l2.ldscore.gz`. Default is `NA`, which will utilize the built-in ld score files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @export
read_ld <- function(ancestry, ld) {
  if (missing(ancestry)) {
    x <- fs::dir_ls(ld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  } else {
    x <- ldscore_files(ancestry, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  }
  return(x)
}
