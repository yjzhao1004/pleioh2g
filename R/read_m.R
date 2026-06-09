#' Read M from either internal or external file - refer to ldscr R package (https://github.com/mglev1n/ldscr)
#'
#' @description
#' `read_m()` Read M from either internal or external file

#' @param ld (character) Path to directory containing ld score files, ending in `*.l2.ldscore.gz`.
#' @return
#' A data frame (tibble) containing SNP counts read from the specified M files.
#'
#' @export
read_m <- function(ld) {

    m <- fs::dir_ls(ld, glob = "*.l2.M_5_50") %>%
      vroom::vroom(col_types = vroom::cols(), col_names = FALSE, delim = "\t")

  return(m)
}
