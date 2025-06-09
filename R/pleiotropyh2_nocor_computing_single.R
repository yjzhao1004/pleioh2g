#' Compute pleioh2g before bias correction for target disease
#'
#' This function is used to compute pleioh2g after bias correction for target disease
#'
#' @param G index of target disease.
#' @param phenotype_path Path to the phenotype name file - a txt file that listed the phenotype names for all GWAS data with header named 'traits'.
#' @param save_path directory where post-corrected pleioh2g results will be stored.
#' @param h2_vector h2g vector for all traits - 1✖number_of_disease matrix - aligned as the order in phenotype file
#' @param h2_vector_mat h2g array from jackknife-block estimates for all traits - blocknum✖number_of_disease matrix - aligned as the order in phenotype file
#' @param Results_full_rg genetic correlation matrix.- number_of_disease✖number_of_disease matrix - aligned as the order in phenotype file
#' @param Results_full_rg_array genetic correlation jackknife-block array.- number_of_disease✖number_of_disease✖blocknum matrix - aligned as the order in phenotype file
#' @import dplyr
#' @import data.table
#' @importFrom stats sd
#' @importFrom utils write.csv
#' @return A `list` containing the following elements:
#'   - `D_T` (character): The value "401.1".
#'   - `percentage_h2pleio` (numeric): The value 0.5.
#'   - `percentage_h2pleio_se` (numeric): The value 0.1.
#'
#' @export
#'
#' @examples
#' G <- 1
#' data(Results_full_rg_15D)
#' data(Results_full_rg_array_15D)
#' data(h2_vector_15D)
#' data(h2_vector_mat_15D)
#' phenotype_path <- system.file("extdata", "phenotype_name_package.txt",
#' package = "pleioh2g")
#' save_path <- system.file("extdata", "save_results_test",package = "pleioh2g")
#' results<-pleiotropyh2_nocor_computing_single(G,phenotype_path,save_path,h2_vector_15D,
#' h2_vector_mat_15D,Results_full_rg_15D,Results_full_rg_array_15D)
#'
pleiotropyh2_nocor_computing_single<-function(G,phenotype_path,save_path,h2_vector,h2_vector_mat,Results_full_rg,Results_full_rg_array){

  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  target_phenotypes <- data.table::fread(phenotype_path)
  target_phenotypes <- as.character(target_phenotypes$traits)

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

  saveRDS(h2pleiotropy_sum_percentage,file = file.path(save_path,paste0(target_diseases_code,'_perpleioh2g.rds')))

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
  Results_collection_category$selected_pleio_pheno<-paste0(colnames(Results_full_rg),collapse = ',')
  Results_collection_category$h2pleiotropy<-h2pleiotropy_sample_sum_allgenome
  Results_collection_category$h2pleiotropy_se<-h2pleiotropy_jkse
  Results_collection_category$percentage_jackknife<-h2pleiotropy_sum_percentage_allgenome
  Results_collection_category$percentage_jackknife_se<-h2pleiotropy_percentage_jkse


  write.csv(Results_collection_category,file = file.path(save_path,paste0(target_diseases_code,'_h2pleio_dresults.csv')),row.names = FALSE)

  return(list(D_T = target_diseases_code, percentage_h2pleio = h2pleiotropy_sum_percentage_allgenome,percentage_h2pleio_se = h2pleiotropy_percentage_jkse))
}

