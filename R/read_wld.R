#' Read wld from either internal or external file
#'
#' @description
#' `read_wld()` Read wld from either internal or external file

#' @param ancestry (character) One of "AFR", "AMR", "CSA", "EAS", "EUR", or "MID", which will utilize the appropriate built-in `ld` and `wld` files from Pan-UK Biobank. If empty or `NULL`, the user must specify paths to `ld` and `wld` files.
#' @param wld (character) Path to directory containing weight files. Default is `NA`, which will utilize the built-in weight files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @export
read_wld <- function(ancestry, wld) {
  if (missing(ancestry)) {
    w <- fs::dir_ls(wld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  } else {
    w <- ldscore_files(ancestry, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  }
  return(w)
}
