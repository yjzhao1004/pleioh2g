#' genomic-block jackknife and compute rg + h2g
#'
#' This function performs genomic-block jackknife and computes rg + h2g.
#'
#' @param n_block number of jackknife blocks.
#' @param hmp3 Directory for hapmap 3 snplist.
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
#' @importFrom stats cov
#' @importFrom stats sd
#' @importFrom utils write.csv
#'
#' @export
#' @examples
#' hmp3 <- system.file("extdata", "w_hm3.snplist",
#' package = "pleioh2g")
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
#' n_block=3
#' Cal_rg_h2g_jk_alltraits(n_block, hmp3, phenotype_path, gwas_munge_dir, save_path, ld_path, wld_path, sample_prev_path, population_prev_path)


Cal_rg_h2g_jk_alltraits <- function(n_block = 200, hmp3, phenotype_path, gwas_munge_dir, save_path, ld_path, wld_path, sample_prev_path = NULL, population_prev_path = NULL) {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  # Load parameters
  ## load phenotype names
  target_phenotypes <- data.table::fread(phenotype_path)
  target_phenotypes <- as.character(target_phenotypes$traits)


  hmp3_snp<-data.table::fread(hmp3)
  blocks <- split(hmp3_snp$SNP, cut(seq_along(hmp3_snp$SNP), breaks = n_block, labels = FALSE))

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

  # create saved gcorr matrix and h2g matrix
  Results_full_h2_array<-matrix(0,nrow = n_block,ncol = length(target_phenotypes))
  rownames(Results_full_h2_array)<-paste0("block_",c(1:n_block))
  colnames(Results_full_h2_array)<-target_phenotypes

  if (!is.null(sample_prev_path) && !is.null(population_prev_path)) {
    Results_full_h2_lia_array<-matrix(0,nrow = n_block,ncol = length(target_phenotypes))
    rownames(Results_full_h2_lia_array)<-paste0("block_",c(1:n_block))
    colnames(Results_full_h2_lia_array)<-target_phenotypes
  }

  Results_full_rg<-matrix(0,nrow = length(target_phenotypes),ncol = length(target_phenotypes))
  rownames(Results_full_rg)<-target_phenotypes
  colnames(Results_full_rg)<-target_phenotypes
  diag(Results_full_rg)<-1
  Results_full_rg_array<-replicate(n_block,Results_full_rg)
  rm(Results_full_rg)


  # start genomic-block jackknife
  for(block in c(1:n_block)){
    print(paste0("jackknife genomic block: ",block))
    block_current <- blocks[[block]]

    # read gwas .sumstat data
    GWAS_list <- lapply(1:length(gwas_munge_dir), function(i) {
      data.table::fread(gwas_munge_dir[i]) %>% dplyr::filter(!is.na(N))
    })

    GWAS_remainblocks <- lapply(GWAS_list, function(gwas) subset(gwas, gwas$SNP %in% setdiff(gwas$SNP, block_current)))

    rg_res <- ldscr::ldsc_rg(
      munged_sumstats = setNames(GWAS_remainblocks, paste0('GWAS_', c(1:length(target_phenotypes)), '_remainblock')),
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

    #save rg matrix
    Results_full_rg<-matrix(0,nrow = length(target_phenotypes),ncol = length(target_phenotypes))
    rownames(Results_full_rg)<-target_phenotypes
    colnames(Results_full_rg)<-target_phenotypes
    lower_tri_indices <- which(lower.tri(Results_full_rg, diag = FALSE), arr.ind = TRUE)
    Results_full_rg[lower_tri_indices]<-rg_res$rg$rg
    Results_full_rg<-Results_full_rg  + t(Results_full_rg)
    diag(Results_full_rg)<-1

    #save gcov matrix
    Results_full_gcov<-rg_res$raw$S
    rownames(Results_full_gcov)<-target_phenotypes
    colnames(Results_full_gcov)<-target_phenotypes


    if (!is.null(sample_prev_path) && !is.null(population_prev_path)) {
      Results_full_h2_lia<-matrix(0,nrow = 1,ncol = length(target_phenotypes))
      colnames(Results_full_h2_lia)<-target_phenotypes
      for(pheno in c(1:length(target_phenotypes))){
        Results_full_h2_lia[1,target_phenotypes[pheno]]<-ldscr::h2_liability(Results_full_h2[pheno], sample_prev$pre[pheno], population_prev$pre[pheno])
      }
    }

    Results_full_h2_array[block,]<-Results_full_h2

    if (!is.null(sample_prev_path) && !is.null(population_prev_path)) {
      Results_full_h2_lia_array[block,]<-Results_full_h2_lia
    }
    Results_full_rg_array[,,block]<-Results_full_rg

  }

  saveRDS(Results_full_h2_array, file.path(save_path,"Results_full_h2_array.rds"))
  if (!is.null(sample_prev_path) && !is.null(population_prev_path)) {
    saveRDS(Results_full_h2_lia_array, file.path(save_path,"Results_full_h2_lia_array.rds"))
  }
  saveRDS(Results_full_rg_array, file.path(save_path,"Results_full_rg_array.rds"))

}
