#' Example munged dataframe - refer to ldscr R package (https://github.com/mglev1n/ldscr)
#'
#' @param dataframe (logical) If `TRUE` (default), return an example munged dataframe. If `FALSE`, return path to the file on disk.
#' @param example (character) "401.1" which have been included as example traits.
#' @return either a [tibble][tibble::tibble-package] containing a munged dataframe, or a path to the file on disk.
#'
#' @export

sumstats_munged_example_input <- function(example, dataframe = TRUE) {
  if (dataframe) {
    vroom::vroom(fs::path(fs::path_package("extdata", paste0(example, ".sumstats.gz"), package = "pleioh2g")), col_types = vroom::cols())
  } else {
    fs::path(fs::path_package("extdata", paste0(example, ".sumstats.gz"), package = "pleioh2g"))
  }
}
