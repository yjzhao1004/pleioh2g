#' Read summary statistics from either internal or external file - refer to ldscr R package (https://github.com/mglev1n/ldscr)
#' @description
#' `read_sumstats()` Read summary statistics from either internal or external file

#' @param munged_sumstats Either a dataframe, or a path to a file containing munged summary statistics. Must contain at least columns named `SNP` (rsid), `A1` (effect allele), `A2` (non-effect allele), `N` (total sample size) and `Z` (Z-score)
#' @param name trait name

#' @export

read_sumstats <- function(munged_sumstats, name) {
  # Check if name is present to enable more informative logging
  if (missing(name)) {
    if (is.character(munged_sumstats)) {
      cli::cli_progress_step("Reading summary statistics from {munged_sumstats}")
      sumstats_df <- vroom::vroom(munged_sumstats, col_types = vroom::cols())
    } else {
      cli::cli_progress_step("Reading summary statistics from dataframe")
      sumstats_df <- munged_sumstats
    }
  } else {
    if (is.character(munged_sumstats)) {
      cli::cli_progress_step("Reading summary statistics for '{name}' from {munged_sumstats}")
      sumstats_df <- vroom::vroom(.x, col_types = vroom::cols())
    } else {
      cli::cli_progress_step("Reading summary statistics for '{name}' from dataframe")
      sumstats_df <- .x
    }
    sumstats_df <- na.omit(sumstats_df)

    return(sumstats_df)
  }
}
