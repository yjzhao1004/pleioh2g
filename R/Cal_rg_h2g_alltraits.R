#' Compute rg + h2g
#'
#' This function is used to compute rg + h2g using LDSC.
#'
#' @param phenotype_path Path to the phenotype name file - a txt file that listed the phenotype names for all GWAS data with header named 'traits'.
#' @param gwas_munge_dir Vector of path for all LDSC-munged GWAS .stat.gz
#' @param save_path Directory where rg + h2g results will be stored.
#' @param ld_path Path to directory containing ld score files.
#' @param wld_path Path to directory containing weight files.
#' @param sample_prev_path Directory for sample prevalence, having two columns, one with header named 'diseases',the other one with header named 'pre'.
#' @param population_prev_path Directory for population prevalence, having two columns, one with header named 'diseases',the other one with header named 'pre'.
#'
#' @import dplyr
#' @import data.table
#' @import ldscr
#' @importFrom dplyr %>%
#' @importFrom stats cov
#' @importFrom stats sd
#' @importFrom utils write.csv
#'
#' @export
#' @examples
#' phenotype_path<-system.file("extdata", "phenotype_name_package_test.txt",
#' package = "pleioh2g")
#' gwas_munge_dir <-  c(system.file("extdata/GWAS_sumstat","401.1.sumstats.gz",
#' package = "pleioh2g"),system.file("extdata/GWAS_sumstat","250.2.sumstats.gz",
#' package = "pleioh2g"),system.file("extdata/GWAS_sumstat","296.22.sumstats.gz",
#' package = "pleioh2g"))
#' save_path <- system.file("extdata", "rg_results_test",package = "pleioh2g")
#' ld_path<-system.file("extdata", "eur_w_ld_chr",package = "pleioh2g")
#' wld_path<-system.file("extdata", "eur_w_ld_chr",package = "pleioh2g")
#' sample_prev_path = system.file("extdata", "samp_prev_test.txt",package = "pleioh2g")
#' population_prev_path = system.file("extdata", "pop_prev_test.txt",package = "pleioh2g")
#' Cal_rg_h2g_alltraits(phenotype_path, gwas_munge_dir, save_path, ld_path, wld_path, sample_prev_path, population_prev_path)

Cal_rg_h2g_alltraits <- function(phenotype_path, gwas_munge_dir, save_path, ld_path, wld_path, sample_prev_path = NULL, population_prev_path = NULL) {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  # Load parameters
  ## load phenotype names
  target_phenotypes <- data.table::fread(phenotype_path)
  target_phenotypes <- as.character(target_phenotypes$traits)


  # if input of prev is null, set default NA
  if (!is.null(sample_prev_path)) {
    sample_prev <-data.table::fread(sample_prev_path)
  }else{
    sample_prev<-NA
  }
  if (!is.null(population_prev_path)) {
    population_prev <-data.table::fread(population_prev_path)
  }else{
    population_prev<-NA
  }


  # read gwas .sumstat data
  GWAS_list <- lapply(1:length(gwas_munge_dir), function(i) {
    data.table::fread(gwas_munge_dir[i]) %>% dplyr::filter(!is.na(N))
  })

  rg_res <- ldscr::ldsc_rg(
    munged_sumstats = setNames(GWAS_list, paste0('GWAS_', c(1:length(target_phenotypes)))),
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

  saveRDS(Results_full_h2, file = file.path(save_path,"Results_full_h2.rds"))
  saveRDS(Results_full_rg, file = file.path(save_path,"Results_full_rg.rds"))
  saveRDS(Results_full_rgz, file = file.path(save_path,"Results_full_rgz.rds"))
  saveRDS(Results_full_gcov, file = file.path(save_path,"Results_full_gcov.rds"))
  saveRDS(Results_full_h2Z, file = file.path(save_path,"Results_full_h2Z.rds"))

  if (!is.null(sample_prev_path) && !is.null(population_prev_path)) {
    Results_full_h2_lia<-matrix(0,nrow = 1,ncol = length(target_phenotypes))
    colnames(Results_full_h2_lia)<-target_phenotypes
    for(pheno in c(1:length(target_phenotypes))){
      Results_full_h2_lia[1,target_phenotypes[pheno]]<-ldscr::h2_liability(Results_full_h2[pheno], sample_prev$pre[pheno], population_prev$pre[pheno])
    }
    saveRDS(Results_full_h2_lia, file.path(save_path,"Results_full_h2_lia.rds"))
  }
}


