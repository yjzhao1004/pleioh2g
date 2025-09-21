#' Read wld from either internal or external file
#'
#' @description
#' `read_wld()` Read wld from either internal or external file
#' @param wld (character) Path to directory containing weight files. Default is `NA`, which will utilize the built-in weight files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @export
read_wld <- function(wld) {

    w <- fs::dir_ls(wld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())

  return(w)
}
