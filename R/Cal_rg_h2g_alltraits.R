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
#'
#' @export
#' @examples
#' phenotype<-c('401.1','250.2','296.22')
#' munged_sumstats = list("401.1" = sumstats_munged_example_input(example = "401.1"), "250.2" = sumstats_munged_example_input(example = "250.2"),"296.22" = sumstats_munged_example_input(example = "296.22"))
#' ld_path<-fs::path(fs::path_package("extdata/eur_w_ld_chr", package = "pleioh2g"))
#' wld_path<-fs::path(fs::path_package("extdata/eur_w_ld_chr", package = "pleioh2g"))
#' sample_prev <- c(0.37,0.1,0.14)
#' population_prev <- c(0.37,0.06,0.14)
#' results<-Cal_rg_h2g_alltraits(phenotype, munged_sumstats, ld_path, wld_path, sample_prev, population_prev)


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
  GWAS_list <- lapply(1:length(munged_sumstats), function(i) {
    munged_sumstats[[i]] %>% dplyr::filter(!is.na(N))
  })

  rg_res <- ldsc_rg(
    munged_sumstats = setNames(GWAS_list, paste0('GWAS_', target_phenotypes)),
    ld = ld_path,
    wld = wld_path,
    n_blocks = 200,
    chisq_max = NA,
    chr_filter = c(1:22)
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
      Results_full_h2_lia[1,target_phenotypes[pheno]]<-ldscr::h2_liability(Results_full_h2[pheno], sample_prev[pheno], population_prev[pheno])
    }
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


