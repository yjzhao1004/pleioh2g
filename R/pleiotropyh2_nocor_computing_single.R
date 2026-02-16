#' Compute pleioh2g before bias correction for target disease
#'
#' This function is used to compute pleioh2g after bias correction for target disease
#'
#' @param G index of target disease.
#' @param phenotype Vector of the phenotype name
#' @param h2_vector h2g vector for all traits  - aligned as the order in phenotype file
#' @param h2_vector_mat h2g array from jackknife-block estimates for all traits - aligned as the order in phenotype file
#' @param Results_full_rg genetic correlation matrix.- aligned as the order in phenotype file
#' @param Results_full_rg_array genetic correlation jackknife-block array.- aligned as the order in phenotype file
#' @import dplyr
#' @import data.table
#' @importFrom stats sd
#' @importFrom utils write.csv
#' @return A `list` containing the following elements:
#'   - `target_disease` (character): The value "401.1".
#'   - `target_disease_h2_est` (numeric): target disease h2g.
#'   - `target_disease_h2_se` (numeric): target disease h2g_se.
#'   - `selected_auxD` (character): auxiliary diseases.
#'   - `h2pleio_uncorr` (numeric): pre-correction pleiotropic heritability estimate.
#'   - `h2pleio_uncorr_se` (numeric): pre-correction pleiotropic heritability jackknife s.e. estimate.
#'   - `percentage_h2pleio_uncorr` (numeric): pre-correction percentage of pleiotropic heritability estimate.
#'   - `percentage_h2pleio_uncorr_se` (numeric): pre-correction percentage of pleiotropic heritability jackknife s.e. estimate.
#'   - `percentage_h2pleio_jackknife_uncorr` (numeric): vector of all pre-correction percentage of pleiotropic heritability jackknife estimates.
#'
#' @export
#'
#' @examples
#' G <- 1
#' data(Results_full_rg_15D)
#' data(Results_full_rg_array_15D)
#' data(h2_vector_15D)
#' data(h2_vector_mat_15D)
#' phenotype<-c("401.1","244.5","318","735.3","411.4",
#' "427.2","454.1","278.1","250.2","550.1","530.11",
#' "296.22","519.8","562.1","763")
#' h2pleiobeforecorr<-pleiotropyh2_nocor_computing_single(G,phenotype,h2_vector_15D,
#' h2_vector_mat_15D,Results_full_rg_15D,Results_full_rg_array_15D)
#'
pleiotropyh2_nocor_computing_single<-function(G,phenotype,h2_vector,h2_vector_mat,Results_full_rg,Results_full_rg_array){

  target_phenotypes <- as.character(phenotype)

  n_block <- dim(Results_full_rg_array)[3]
  n_trait <- dim(Results_full_rg)[1]

  target_diseases_code<-target_phenotypes[G]

  Results_collection_category=list()

  #computing h2pleio for jk blocks
  h2pleiotropy_sample_sum<- vector()
  targetnum<-which(colnames(Results_full_rg_array)==target_diseases_code)
  h2total_target_allgenome<-h2_vector[targetnum]
  h2total_target<-h2_vector_mat[,targetnum]

  for (t in c(1:n_block)){
    h2pleiotropy_sample_sum[t] <- Cal_cor_pleiotropic_h2_single(Results_full_rg_array[,,t],h2total_target[t],targetnum)
  }
  h2pleiotropy_sum_percentage<-h2pleiotropy_sample_sum/h2total_target

  #saveRDS(h2pleiotropy_sum_percentage,file = file.path(save_path,paste0(target_diseases_code,'_perpleioh2g.rds')))

  #computing pleioh2g (point est)
  h2pleiotropy_sample_sum_allgenome<-Cal_cor_pleiotropic_h2_single(Results_full_rg,h2total_target_allgenome,targetnum)
  h2pleiotropy_sum_percentage_allgenome<-h2pleiotropy_sample_sum_allgenome/h2total_target_allgenome

  #computing pleioh2g (jk se)
  pseudovalue_pleioh2g <- matrix(data = NA, nrow = n_block, ncol = 1)
  pseudovalue_h2g <- matrix(data = NA, nrow = n_block, ncol = 1)
  pseudovalue_per_pleioh2g <- matrix(data = NA, nrow = n_block, ncol = 1)
  for (b in c(1:n_block)){
    pseudovalue_pleioh2g[b,]<-n_block*h2pleiotropy_sample_sum_allgenome-(n_block-1)*h2pleiotropy_sample_sum[b]
    pseudovalue_h2g[b,]<-n_block*h2total_target_allgenome-(n_block-1)*h2total_target[b]
    pseudovalue_per_pleioh2g[b,]<-n_block*h2pleiotropy_sum_percentage_allgenome-(n_block-1)*h2pleiotropy_sum_percentage[b]
  }
  h2_target_jkse<-sd(pseudovalue_h2g)/sqrt(n_block)
  pseudovalue_pleioh2g_mean<-mean(pseudovalue_pleioh2g)
  pseudovalue_h2g_mean<-mean(pseudovalue_h2g)
  pseudovalue_per_pleioh2g_mean<-mean(pseudovalue_per_pleioh2g)
  h2pleiotropy_jkse<-sd(pseudovalue_pleioh2g)/sqrt(n_block)
  h2pleiotropy_percentage_jkse<-sd(pseudovalue_per_pleioh2g)/sqrt(n_block)

  #collect results
  Results_collection_category$target_disease<-target_diseases_code
  Results_collection_category$target_disease_h2_est<-h2total_target_allgenome
  Results_collection_category$target_disease_h2_se<-h2_target_jkse
  Results_collection_category$selected_auxD<-paste0(colnames(Results_full_rg),collapse = ',')
  Results_collection_category$h2pleio_uncorr<-h2pleiotropy_sample_sum_allgenome
  Results_collection_category$h2pleio_uncorr_se<-h2pleiotropy_jkse
  Results_collection_category$percentage_h2pleio_uncorr<-h2pleiotropy_sum_percentage_allgenome
  Results_collection_category$percentage_h2pleio_uncorr_se<-h2pleiotropy_percentage_jkse
  Results_collection_category$percentage_h2pleio_jackknife_uncorr<-h2pleiotropy_sum_percentage

  return(Results_collection_category)
}

