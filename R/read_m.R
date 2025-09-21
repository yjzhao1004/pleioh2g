#' Read M from either internal or external file
#'
#' @description
#' `read_m()` Read M from either internal or external file

#' @param ld (character) Path to directory containing ld score files, ending in `*.l2.ldscore.gz`. Default is `NA`, which will utilize the built-in ld score files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @export
read_m <- function(ld) {

    m <- fs::dir_ls(ld, glob = "*.l2.M_5_50") %>%
      vroom::vroom(col_types = vroom::cols(), col_names = FALSE, delim = "\t")

  return(m)
}
