#' genomic-block jackknife and compute rg + h2g
#'
#' This function performs genomic-block jackknife and computes rg + h2g.
#'
#' @param n_block number of jackknife blocks.
#' @param hmp3 Directory for hapmap 3 snplist.
#' @param phenotype Vector of the phenotype name
#' @param munged_sumstats All LDSC-munged GWAS .stat.gz
#' @param ld_path Path to directory containing ld score files.
#' @param wld_path Path to directory containing weight files.
#' @param sample_prev Vector of sample prevalence, in the same order of input GWAS summary statistics.
#' @param population_prev Vector of population prevalence, in the same order of input GWAS summary statistics.
#'
#' @import dplyr
#' @import data.table
#' @importFrom stats cov
#' @importFrom stats sd
#' @importFrom utils write.csv
#' @return
#' A named list containing block jackknife estimates of SNP-heritability and genetic
#' correlation across all input phenotypes. The list includes the following elements:
#' \itemize{
#'   \item \code{h2array}: A matrix of per-block SNP-heritability estimates on the
#'         observed scale. Rows correspond to jackknife blocks, and columns correspond
#'         to input phenotypes.
#'   \item \code{liah2array}: A matrix of per-block SNP-heritability estimates on the
#'         liability scale, with the same row and column structure as \code{h2array}.
#'   \item \code{rgarray}: A three-dimensional array of pairwise genetic correlation
#'         estimates. The first two dimensions represent phenotype pairs
#'         (rows and columns), and the third dimension indexes the jackknife blocks.
#'   \item \code{gcovarray}: A three-dimensional array of pairwise genetic covariance
#'         estimates, aligned in structure with \code{rgarray}.
#' }
#'
#' Each element provides per-block estimates that can be used to compute
#' standard errors or confidence intervals via the block jackknife method.
#'
#' @export
#'




Cal_rg_h2g_jk_alltraits <- function(n_block = 200, hmp3, phenotype, munged_sumstats, ld_path, wld_path, sample_prev = NULL, population_prev = NULL) {
  ## load phenotype names
  target_phenotypes <- as.character(phenotype)

  hmp3_snp<-data.table::fread(hmp3)
  blocks <- split(hmp3_snp$SNP, cut(seq_along(hmp3_snp$SNP), breaks = n_block, labels = FALSE))

  # if input of prev is null, set default NA
  if (is.null(sample_prev)) {
    sample_prev<-NA
  }
  if (is.null(population_prev)) {
    population_prev<-NA
  }

  # create saved gcorr matrix and h2g matrix
  Results_full_h2_array<-matrix(0,nrow = n_block,ncol = length(target_phenotypes))
  rownames(Results_full_h2_array)<-paste0("block_",c(1:n_block))
  colnames(Results_full_h2_array)<-target_phenotypes

  if (!is.null(sample_prev) && !is.null(population_prev)) {
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

  Results_full_gcov<-matrix(0,nrow = length(target_phenotypes),ncol = length(target_phenotypes))
  rownames(Results_full_gcov)<-target_phenotypes
  colnames(Results_full_gcov)<-target_phenotypes
  Results_full_gcov_array<-replicate(n_block,Results_full_gcov)
  rm(Results_full_gcov)


  # start genomic-block jackknife
  for(block in c(1:n_block)){
    message(paste0("jackknife genomic block: ",block))
    block_current <- blocks[[block]]

    # read gwas .sumstat data
    GWAS_list <- lapply(munged_sumstats, function(df) {
      dplyr::filter(df, !is.na(.data$N))   # 用 .data$N 避免全局变量警告
    })

    GWAS_remainblocks <- lapply(GWAS_list, function(gwas) {
      subset(gwas, SNP %in% setdiff(gwas$SNP, block_current))
    })

    rg_res <- ldsc_rg(
      munged_sumstats = stats::setNames(
        GWAS_remainblocks,
        paste0("GWAS_", target_phenotypes, "_remainblock")
      ),
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


    if (!is.null(sample_prev) && !is.null(population_prev)) {
      Results_full_h2_lia<-matrix(0,nrow = 1,ncol = length(target_phenotypes))
      colnames(Results_full_h2_lia)<-target_phenotypes
      for(pheno in c(1:length(target_phenotypes))){
        Results_full_h2_lia[1,target_phenotypes[pheno]]<-h2_liability(Results_full_h2[pheno], sample_prev[pheno], population_prev[pheno])
      }
      message("Estimated liability scale heritability..done.")
    }

    Results_full_h2_array[block,]<-Results_full_h2

    if (!is.null(sample_prev) && !is.null(population_prev)) {
      Results_full_h2_lia_array[block,]<-Results_full_h2_lia
    }
    Results_full_rg_array[,,block]<-Results_full_rg
    Results_full_gcov_array[,,block]<-Results_full_gcov
  }

  output<-list(
    h2array = Results_full_h2_array,
    rgarray = Results_full_rg_array,
    gcovarray = Results_full_gcov_array
  )

  if (!is.null(sample_prev) && !is.null(population_prev)) {
    output<-list(
      h2array = Results_full_h2_array,
      liah2array = Results_full_h2_lia_array,
      rgarray = Results_full_rg_array,
      gcovarray = Results_full_gcov_array
    )

  }
  return(output)
}
