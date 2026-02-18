#' Compute rg + h2g
#'
#' This function is used to compute rg + h2g using LDSC.
#'
#' @param phenotype Vector of the phenotype name
#' @param munged_sumstats All LDSC-munged GWAS .stat.gz
#' @param ld_path Path to directory containing ld score files.
#' @param wld_path Path to directory containing weight files.
#' @param sample_prev Vector of sample prevalence, in the same order of input GWAS summary statistics.
#' @param population_prev Vector of population prevalence, in the same order of input GWAS summary statistics.
#'
#' @import dplyr
#' @import data.table
#' @importFrom dplyr %>%
#' @importFrom stats cov
#' @importFrom stats sd
#' @importFrom utils write.csv
#' @return
#' A named list containing LDSC-based heritability and genetic correlation estimates
#' across all input phenotypes. The list includes the following elements:
#' \itemize{
#'   \item \code{h2}: Matrix of SNP-heritability estimates on the observed scale
#'         (rows = 1, columns = input phenotypes).
#'   \item \code{h2Z}: Matrix of corresponding heritability Z-scores.
#'   \item \code{liah2}: Matrix of heritability estimates on the liability scale.
#'   \item \code{rg}: Symmetric matrix of pairwise genetic correlations between traits.
#'   \item \code{rgz}: Matrix of Z-scores for the genetic correlation estimates.
#'   \item \code{gcov}: Symmetric matrix of genetic covariances between traits.
#' }
#'
#' Each element corresponds to one LDSC-derived summary statistic, with trait names
#' used as both row and column names.
#'
#'
#' @export



Cal_rg_h2g_alltraits <- function(phenotype, munged_sumstats, ld_path, wld_path, sample_prev = NULL, population_prev = NULL) {

  ## load phenotype names
  target_phenotypes <- as.character(phenotype)

  # if input of prev is null, set default NA
  if (is.null(sample_prev)) {
    sample_prev<-NA
  }
  if (is.null(population_prev)) {
    population_prev<-NA
  }

  # read gwas .sumstat data
  GWAS_list <- lapply(munged_sumstats, function(df) {
    dplyr::filter(df, !is.na(.data$N))      # 使用 .data$N 避免全局变量检查
  })

  rg_res <- ldsc_rg(
    munged_sumstats = stats::setNames(GWAS_list, paste0("GWAS_", target_phenotypes)),
    ld = ld_path,
    wld = wld_path,
    n_blocks = 200,
    chisq_max = NA,
    chr_filter = 1:22
  )
  # create saved gcorr matrix and h2g matrix

  #save observed h2
  Results_full_h2<-matrix(0,nrow = 1,ncol = length(target_phenotypes))
  colnames(Results_full_h2)<-target_phenotypes
  Results_full_h2[1,]<-rg_res$h2$h2_observed

  Results_full_h2Z<-matrix(0,nrow = 1,ncol = length(target_phenotypes))
  colnames(Results_full_h2Z)<-target_phenotypes
  Results_full_h2Z[1,]<-rg_res$h2$h2_Z

  #save rg matrix
  Results_full_rg<-matrix(0,nrow = length(target_phenotypes),ncol = length(target_phenotypes))
  rownames(Results_full_rg)<-target_phenotypes
  colnames(Results_full_rg)<-target_phenotypes
  lower_tri_indices <- which(lower.tri(Results_full_rg, diag = FALSE), arr.ind = TRUE)
  Results_full_rg[lower_tri_indices]<-rg_res$rg$rg
  Results_full_rg<-Results_full_rg  + t(Results_full_rg)
  diag(Results_full_rg)<-1

  #save rgz matrix
  Results_full_rgz<-matrix(0,nrow = length(target_phenotypes),ncol = length(target_phenotypes))
  rownames(Results_full_rgz)<-target_phenotypes
  colnames(Results_full_rgz)<-target_phenotypes
  lower_tri_indices <- which(lower.tri(Results_full_rgz, diag = FALSE), arr.ind = TRUE)
  Results_full_rgz[lower_tri_indices]<-rg_res$rg$rg/rg_res$rg$rg_se
  Results_full_rgz<-Results_full_rgz  + t(Results_full_rgz)
  diag(Results_full_rgz)<-NA

  #save gcov matrix
  Results_full_gcov<-rg_res$raw$S
  rownames(Results_full_gcov)<-target_phenotypes
  colnames(Results_full_gcov)<-target_phenotypes

  output<-list(
    h2 = Results_full_h2,
    h2Z = Results_full_h2Z,
    rg = Results_full_rg,
    rgz = Results_full_rgz,
    gcov = Results_full_gcov
  )

  if (!is.null(sample_prev) && !is.null(population_prev)) {

    Results_full_h2_lia<-matrix(0,nrow = 1,ncol = length(target_phenotypes))
    colnames(Results_full_h2_lia)<-target_phenotypes
    for(pheno in c(1:length(target_phenotypes))){
      Results_full_h2_lia[1,target_phenotypes[pheno]]<-h2_liability(Results_full_h2[pheno], sample_prev[pheno], population_prev[pheno])
    }
    message("Estimated liability scale heritability..done.")
    output<-list(
      h2 = Results_full_h2,
      h2Z = Results_full_h2Z,
      liah2 = Results_full_h2_lia,
      rg = Results_full_rg,
      rgz = Results_full_rgz,
      gcov = Results_full_gcov
    )
  }
  return(output)
}


