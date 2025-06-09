#' Compute pleioh2g using all single auxiliary categories
#'
#' @param G index of target disease.
#' @param singleCat single auxiliary categort name - the same as category information in phenotype file
#' @param phenotype_path Path to the phenotype name and h2z score file - 1st col: 'traits'; 2nd col: 'h2_z'
#' @param save_path directory where post-corrected pleioh2g results will be stored.
#' @param h2_vector h2g vector for all traits - 1✖number_of_disease matrix - aligned as the order in phenotype file
#' @param h2_vector_mat h2g array from jackknife-block estimates for all traits - blocknum✖number_of_disease matrix - aligned as the order in phenotype file
#' @param Results_full_rg genetic correlation matrix.- number_of_disease✖number_of_disease matrix - aligned as the order in phenotype file
#' @param Results_full_rg_array genetic correlation jackknife-block array.- number_of_disease✖number_of_disease✖blocknum matrix - aligned as the order in phenotype file
#' @param sample_rep sampling times in bias correction
#' @param tolerance tolerance for dichotomic search
#' @param seed random seed
#' @import dplyr
#' @import data.table
#' @importFrom stats sd
#' @importFrom utils write.csv
#' @export
#'
#' @examples
#' G = 3
#' data("Results_full_rg_45D")
#' data("Results_full_rg_array_45D")
#' data("h2_vector_mat_45D")
#' data("h2_vector_45D")
#' data("Rg_mat_z_45D")
#' phenotype_path <- system.file("extdata", "Disease_45d_h2Z_Category.csv",
#' package = "pleioh2g")
#' save_path <- system.file("extdata", "save_results_test",package = "pleioh2g")
#' sample_rep <- 200
#' tolerance <- 1e-6
#' seed <- 123
#' singleCat<-"mental disorders"
#' pleiotropyh2_singlecategory_single(G,phenotype_path,singleCat,save_path,h2_vector_45D,h2_vector_mat_45D,Results_full_rg_45D,
#' Results_full_rg_array_45D,Rg_mat_z_45D,sample_rep,tolerance,seed)

pleiotropyh2_singlecategory_single<-function(G,phenotype_path,singleCat,save_path,h2_vector,h2_vector_mat,Results_full_rg,Results_full_rg_array,Rg_mat_z,sample_rep,tolerance,seed){
  target_phenotype<-data.table::fread(phenotype_path)
  unique_category<-unique(target_phenotype$Category)
  currentT<-target_phenotype$traits[G]
  currentT_C<-target_phenotype$Category[G]

  currentcategory<-singleCat
  print(paste0("Process: ",currentcategory))
  select_disease<-target_phenotype$traits[target_phenotype$Category%in%currentcategory]

  if(currentcategory==currentT_C){
    select_disease<-select_disease

  }else{
    select_disease<-c(currentT,select_disease)

  }

  save_path_update<-paste0(save_path,'/single_',currentcategory)
  dir.create(save_path_update)

  phenotypename_update<-data.frame(
    traits = select_disease,
    h2_z = target_phenotype$h2_z[match(select_disease,target_phenotype$traits)]
  )
  phenotype_path_update<-paste0(save_path_update,'/',currentT,'_initialpheno.txt')

  data.table::fwrite(phenotypename_update,file = phenotype_path_update,quote = F,sep = '\t')

  h2_vector_update<-t(as.matrix(h2_vector[1,phenotypename_update$traits]))
  h2_vector_mat_update<-h2_vector_mat[,phenotypename_update$traits]
  Results_full_rg_update<-Results_full_rg[phenotypename_update$traits,phenotypename_update$traits]
  Results_full_rg_array_update<-Results_full_rg_array[phenotypename_update$traits,phenotypename_update$traits,]
  Rg_mat_z_update<-Rg_mat_z[phenotypename_update$traits, phenotypename_update$traits]

  target_num = which(phenotypename_update$traits==currentT)
  pruning_pleioh2g_corr_single_rgzscore_cutatall_more(target_num,phenotype_path_update,save_path_update,h2_vector_update,h2_vector_mat_update,Results_full_rg_update,Results_full_rg_array_update,Rg_mat_z_update,sample_rep,tolerance,seed)

}
