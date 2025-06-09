
#' Transform GWAS summary statistics to pre LDSC .stat.gz data for all phenotypes
#'
#' This function transform multiple GWAS summary statistics into .sumstats.gz form by creating necessary directories and running the munge function on each file.
#'
#' @param hmp3 Path to the hapmap3 SNP list file.
#' @param gwas_dir vector of GWAS summary statistics path
#' @param Nsamp vector of sample size - the same order as GWAS in gwas_dir
#' @param trait_names vector of phenotype name for prefix name of.sumstats.gz data - the same order as GWAS in gwas_dir
#' @param output_dir Directory where  LDSC munged format GWAS .sumstats.gz will be stored.
#'
#' @import dplyr
#' @import data.table
#' @import ldscr
#' @import GenomicSEM
#'
#' @export
#' @examples
#' hmp3 <- system.file("extdata", "w_hm3.snplist",
#' package = "pleioh2g")
#' gwas_dir <- c(system.file("extdata/GWAS_example","401.1_gwas.txt.gz",
#' package = "pleioh2g"),system.file("extdata/GWAS_example","250.2_gwas.txt.gz",
#' package = "pleioh2g"),system.file("extdata/GWAS_example","296.22_gwas.txt.gz",
#' package = "pleioh2g"))
#' Nsamp <- c(157206,157206,157206)
#' output_dir <- system.file("extdata/GWAS_sumstat",
#' package = "pleioh2g")
#' trait_names<-c("401.1","250.2","296.22")
#' munge_gwas_allphenotype(hmp3, gwas_dir, Nsamp, trait_names,output_dir)


munge_gwas_allphenotype <- function(hmp3, gwas_dir, Nsamp, trait_names,output_dir) {

  # Check if target_phenotypes and Nsamp have the same length
  if (length(gwas_dir) != length(Nsamp)){
    stop("Error: target_phenotypes and Nsamp must have the same length.")
  }

  # Check if output directory exists, if not, create it
  if (!dir.exists(output_dir)){
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Process each phenotype
  for (num in c(1:length(gwas_dir))) {
    setwd(paste0(output_dir))
    GenomicSEM::munge(gwas_dir[num], hmp3, trait.names=trait_names[num], Nsamp[num])
  }

}


