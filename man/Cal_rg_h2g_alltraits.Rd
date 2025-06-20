% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cal_rg_h2g_alltraits.R
\name{Cal_rg_h2g_alltraits}
\alias{Cal_rg_h2g_alltraits}
\title{Compute rg + h2g}
\usage{
Cal_rg_h2g_alltraits(
  phenotype_path,
  gwas_munge_dir,
  save_path,
  ld_path,
  wld_path,
  sample_prev_path = NULL,
  population_prev_path = NULL
)
}
\arguments{
\item{phenotype_path}{Path to the phenotype name file - a txt file that listed the phenotype names for all GWAS data with header named 'traits'.}

\item{gwas_munge_dir}{Vector of path for all LDSC-munged GWAS .stat.gz}

\item{save_path}{Directory where rg + h2g results will be stored.}

\item{ld_path}{Path to directory containing ld score files.}

\item{wld_path}{Path to directory containing weight files.}

\item{sample_prev_path}{Directory for sample prevalence, having two columns, one with header named 'diseases',the other one with header named 'pre'.}

\item{population_prev_path}{Directory for population prevalence, having two columns, one with header named 'diseases',the other one with header named 'pre'.}
}
\description{
This function is used to compute rg + h2g using LDSC.
}
\examples{
phenotype_path<-system.file("extdata", "phenotype_name_package_test.txt",
package = "pleioh2g")
gwas_munge_dir <-  c(system.file("extdata/GWAS_sumstat","401.1.sumstats.gz",
package = "pleioh2g"),system.file("extdata/GWAS_sumstat","250.2.sumstats.gz",
package = "pleioh2g"),system.file("extdata/GWAS_sumstat","296.22.sumstats.gz",
package = "pleioh2g"))
save_path <- system.file("extdata", "rg_results_test",package = "pleioh2g")
ld_path<-system.file("extdata", "eur_w_ld_chr",package = "pleioh2g")
wld_path<-system.file("extdata", "eur_w_ld_chr",package = "pleioh2g")
sample_prev_path = system.file("extdata", "samp_prev_test.txt",package = "pleioh2g")
population_prev_path = system.file("extdata", "pop_prev_test.txt",package = "pleioh2g")
Cal_rg_h2g_alltraits(phenotype_path, gwas_munge_dir, save_path, ld_path, wld_path, sample_prev_path, population_prev_path)
}
