#' Merging summary statistics with LD-score files
#'
#' @description
#' `merge_sumstats()` Merging summary statistics with LD-score files

#' @param sumstats_df dataframe of sumstat
#' @param w wld score
#' @param x ld score
#' @param chr_filter (numeric vector) Chromosomes to include in analysis. Separating even/odd chromosomes may be useful for exploratory/confirmatory factor analysis.
#'
#' @export

merge_sumstats <- function(sumstats_df, w, x, chr_filter) {
  merged <- sumstats_df %>%
    arrow::as_arrow_table() %>%
    # dtplyr::lazy_dt() %>%
    dplyr::select(SNP, N, Z, A1) %>%
    dplyr::inner_join((w[, c("SNP", "wLD")] %>% arrow::as_arrow_table()), by = c("SNP")) %>%
    dplyr::inner_join((x %>% arrow::as_arrow_table()), by = c("SNP")) %>%
    # dplyr::arrange(CHR, BP) %>%
    dplyr::filter(CHR %in% chr_filter) %>%
    na.omit() %>%
    unique() %>%
    dplyr::collect() %>%
    tibble::as_tibble()

  return(merged)
}
