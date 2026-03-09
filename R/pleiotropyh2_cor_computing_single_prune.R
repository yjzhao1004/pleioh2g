#' Compute pleioh2g after bias correction for target disease
#'
#' This function is used to compute pleioh2g after bias correction for target disease
#'
#' @param G index of target disease.
#' @param phenotype Vector of the phenotype name
#' @param h2_vector h2g vector for all traits - aligned as the order in phenotype file
#' @param h2_vector_mat h2g array from jackknife-block estimates for all traits - aligned as the order in phenotype file
#' @param Results_full_rg genetic correlation matrix. - aligned as the order in phenotype file
#' @param Results_full_rg_array genetic correlation jackknife-block array. - aligned as the order in phenotype file
#' @param sample_rep sampling times in bias correction



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
#'   - `percentage_h2pleio_uncorr_jackknife` (numeric): vector of all pre-correction percentage of pleiotropic heritability jackknife estimates.
#'   - `h2pleio_corr` (numeric): post-correction pleiotropic heritability estimate.
#'   - `h2pleio_corr_se` (numeric): post-correction pleiotropic heritability estimate s.e..
#'   - `percentage_h2pleio_corr` (numeric): post-correction percentage of pleiotropic heritability estimate.
#'   - `percentage_h2pleio_corr_se` (numeric): post-correction percentage of pleiotropic heritability jackknife s.e. estimate.
#'   - `percentage_h2pleio_corr_Z` (numeric): post-correction percentage of pleiotropic heritability estimate Z score.
#'   - `corrected_weight` (numeric): corrected weight in bias correction.
#'
#'
#' @export
#'
#' @examples
#' \donttest{
#' G <- 1
#' data(Results_full_rg)
#' data(Results_full_rg_array)
#' data(h2_vector)
#' data(h2_vector_mat)
#' Results_full_rg<-Results_full_rg[1:15,1:15]
#' Results_full_rg_array<-Results_full_rg_array[1:15,1:15,]
#' h2_vector<-t(as.matrix(h2_vector[1,1:15]))
#' h2_vector_mat<-h2_vector_mat[,1:15]
#' phenotype<-c("401.1","244.5","318","735.3","411.4",
#' "427.2","454.1","278.1","250.2","550.1","530.11",
#' "296.22","519.8","562.1","763")
#' sample_rep<-20
#' post_corrrresults_prune<-pleiotropyh2_cor_computing_single_prune(G,phenotype,h2_vector,
#' h2_vector_mat,Results_full_rg,Results_full_rg_array, sample_rep)
#'}
pleiotropyh2_cor_computing_single_prune<-function(G,phenotype,h2_vector,h2_vector_mat,Results_full_rg,Results_full_rg_array,sample_rep){

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
  Results_collection_category$percentage_h2pleio_uncorr_jackknife<-h2pleiotropy_sum_percentage

  ####correct
  true_pleio_cor_test <- Cal_cor_test_single(Results_full_rg,targetnum)

  if(h2pleiotropy_sum_percentage_allgenome==true_pleio_cor_test){
    print("After check: truth = estimation")
  }else{
    print("error!")
  }

  corg2dg1_mat<-matrix(0,nrow = n_trait,ncol = n_trait-1)
  lower_bound <- matrix(0,nrow = n_trait,ncol = n_trait-1)
  upper_bound <- matrix(0,nrow = n_trait,ncol = n_trait-1)

  for(plei_h2_idx in 1:n_trait){
    g2dg1 <- Results_full_rg[plei_h2_idx, ][-plei_h2_idx]
    upper_bound[plei_h2_idx,] <- g2dg1
    corg2dg1_mat[plei_h2_idx,] <- g2dg1
  }

  upper_bound_single<-upper_bound[targetnum,]
  lower_bound_single<-lower_bound[targetnum,]
  a_diseases<-corg2dg1_mat[targetnum,]
  corg2dg1<-corg2dg1_mat[targetnum,]


  sample_pleio_cor <- matrix(NA, nrow = sample_rep, ncol = 1)
  set.seed(123)
  print("Start computing initial value...")
  for(sample_id in 1:sample_rep){
    # print(sample_id)
    sample_pleio_cor[sample_id,1]  <- generate_proposal_sample_changea_cor(Results_full_rg,Results_full_rg_array,targetnum, 1)
    if(is.na(sample_pleio_cor[sample_id,1])){
      warning("Warning: Computation took more than 200 iterations to sample a non-singular matrix. Stopping execution. Require further pruning..")
      break
    }
  }

  if (is.na(sample_pleio_cor[sample_id, 1])) {
    stop("Error: sample_pleio_cor[sample_id, 1] is NA. Function execution stopped.")
  }

  mean_est_a_single <- apply(sample_pleio_cor, 2, function(x) mean(x, na.rm = T))
  sd_est_a_single <- apply(sample_pleio_cor, 2, function(x) sd(x, na.rm = T))

  print("Finish computing initial value...")

  t=1

  mean_est_a_single_draw<-vector()
  mean_est_a_single_draw<-c(mean_est_a_single_draw,mean_est_a_single)
  sd_est_a_single_draw<-vector()
  sd_est_a_single_draw<-c(sd_est_a_single_draw,sd_est_a_single)
  ratio_a_vector<-vector()
  ratio_a_vector<-c(ratio_a_vector,1)
  print(paste0("mean estimate across MC samples: ",as.numeric(mean_est_a_single)))
  print(paste0("sd across MC samples: ",as.numeric(sd_est_a_single)))
  print(paste0("observed true estimate: ",as.numeric(true_pleio_cor_test)))
  print(paste0("observed true estimate jackknife s.e.: ",as.numeric(h2pleiotropy_percentage_jkse)))
  print(paste0("corrected weight: ",1))
  print(paste0("difference: ",as.numeric(abs(mean_est_a_single-true_pleio_cor_test)/true_pleio_cor_test)))
  print("Start binary search...")

  while (max(abs(upper_bound_single-lower_bound_single)) > 1e-6) {
    print(paste0("search_round:",t))

    if (mean_est_a_single < true_pleio_cor_test) {

      lower_bound_single <- a_diseases
      a_diseases <- (lower_bound_single + upper_bound_single)/2
      ratio_a<-(a_diseases/corg2dg1)[1]

      sample_pleio_a_test <- matrix(NA, nrow = sample_rep, ncol = 1)
      set.seed(123)
      for(sample_id in 1:sample_rep){
        sample_pleio_a_test[sample_id,1]  <- generate_proposal_sample_changea_cor(Results_full_rg,Results_full_rg_array,targetnum, ratio_a)
      }
      mean_est_a_single <- apply(sample_pleio_a_test, 2, function(x) mean(x, na.rm = T))
      mean_est_a_single_draw<-c(mean_est_a_single_draw,mean_est_a_single)
      sd_est_a_single <- apply(sample_pleio_a_test, 2, function(x) sd(x, na.rm = T))
      sd_est_a_single_draw<-c(sd_est_a_single_draw,sd_est_a_single)
      ratio_a_vector<-c(ratio_a_vector,ratio_a)
    }else{
      upper_bound_single <- a_diseases
      a_diseases <- (lower_bound_single + upper_bound_single)/2
      ratio_a<-(a_diseases/corg2dg1)[1]

      sample_pleio_a_test <- matrix(NA, nrow = sample_rep, ncol = 1)
      set.seed(123)
      for(sample_id in 1:sample_rep){
        sample_pleio_a_test[sample_id,1]  <- generate_proposal_sample_changea_cor(Results_full_rg,Results_full_rg_array,targetnum, ratio_a)
      }
      mean_est_a_single <- apply(sample_pleio_a_test, 2, function(x) mean(x, na.rm = T))
      mean_est_a_single_draw<-c(mean_est_a_single_draw,mean_est_a_single)
      sd_est_a_single <- apply(sample_pleio_a_test, 2, function(x) sd(x, na.rm = T))
      sd_est_a_single_draw<-c(sd_est_a_single_draw,sd_est_a_single)
      ratio_a_vector<-c(ratio_a_vector,ratio_a)
    }

    print(paste0("mean estimate across MC samples: ",as.numeric(mean_est_a_single)))
    print(paste0("sd across MC samples: ",as.numeric(sd_est_a_single)))
    print(paste0("observed true estimate: ",as.numeric(true_pleio_cor_test)))
    print(paste0("observed true estimate jackknife s.e.: ",as.numeric(h2pleiotropy_percentage_jkse)))
    print(paste0("corrected weight: ",as.numeric(ratio_a)))
    print(paste0("difference: ",as.numeric(abs(mean_est_a_single-true_pleio_cor_test)/true_pleio_cor_test)))


    if(abs(mean_est_a_single-true_pleio_cor_test)/true_pleio_cor_test<0.05){
      break
    }

    t=t+1
  }

  corrected_weightd<-ratio_a

  true_pleio_h2_correct<-Cal_cor_pleiotropic_h2_corrected_single(Results_full_rg,h2total_target_allgenome,corrected_weightd,targetnum)
  true_pleio_h2_correct_percentage<-true_pleio_h2_correct/h2total_target_allgenome


  Results_collection_category$h2pleio_corr<-true_pleio_h2_correct
  post_se<-sd_est_a_single/sd_est_a_single_draw[1]*h2pleiotropy_percentage_jkse
  Results_collection_category$h2pleio_corr_se<-post_se*h2total_target_allgenome

  Results_collection_category$percentage_h2pleio_corr<-true_pleio_h2_correct_percentage
  Results_collection_category$percentage_h2pleio_corr_se<-sd_est_a_single/sd_est_a_single_draw[1]*h2pleiotropy_percentage_jkse
  Results_collection_category$percentage_h2pleio_corr_Z<-Results_collection_category$percentage_h2pleio_corr/Results_collection_category$percentage_h2pleio_corr_se

  Results_collection_category$corrected_weight<-corrected_weightd

  return(Results_collection_category)
}
