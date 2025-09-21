#' Read ld from either internal or external file
#'
#' @description
#' `read_ld()` Read ld from either internal or external file.
#' @param ld (character) Path to directory containing ld score files, ending in `*.l2.ldscore.gz`. Default is `NA`, which will utilize the built-in ld score files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @export
read_ld <- function(ld) {

    x <- fs::dir_ls(ld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())

  return(x)
}
