#' Compute pleioh2g excluding target disease category
#'
#' @param G index of target disease.
#' @param phenotype_path Path to the phenotype name and h2z score file - 1st col: 'traits'; 2nd col: 'h2_z'
#' @param allauxD_T_results_path path to the result using all auxiliary disease categories but removing Target category - the same target disease as G
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
#' save_path <- system.file("extdata", "save_results_test/D_T_others",package = "pleioh2g")
#' sample_rep <- 10
#' tolerance <- 1e-6
#' seed <- 123
#' allauxD_T_results_path<-system.file("extdata/save_results_test/D_T", "318_h2pleio_dresults.csv",
#' package = "pleioh2g")
#' pleiotropyh2_D_T_othersall(G,phenotype_path,save_path,allauxD_T_results_path,h2_vector_45D,h2_vector_mat_45D,Results_full_rg_45D,
#' Results_full_rg_array_45D,Rg_mat_z_45D,sample_rep,tolerance,seed)

pleiotropyh2_D_T_othersall<-function(G,phenotype_path,save_path,allauxD_T_results_path,h2_vector,h2_vector_mat,Results_full_rg,Results_full_rg_array,Rg_mat_z,sample_rep,tolerance,seed){
  target_phenotype<-data.table::fread(phenotype_path)
  unique_category<-unique(target_phenotype$Category)
  currentT<-target_phenotype$traits[G]
  currentT_C<-target_phenotype$Category[G]

  unique_category_left<-setdiff(unique_category,currentT_C)

  file_DG<-data.table::fread(allauxD_T_results_path)
  pruneddisease<- unlist(strsplit(file_DG$selected_pleio_pheno, ','))
  pruneddisease_category<-subset(target_phenotype,traits%in%pruneddisease)

  current_exT <- setdiff(pruneddisease_category$traits, pruneddisease_category$traits[which(pruneddisease_category$Category == currentT_C)])

  for(C in c(1:length(unique_category_left))){
    current_reC<-unique_category_left[C]
    save_path_update<-paste0(save_path,'/de_',unique_category_left[C])
    dir.create(save_path_update)

    current_exT_exC<-setdiff(current_exT,pruneddisease_category$traits[which(pruneddisease_category$Category==current_reC)])

    phenotypename_update<-data.frame(
      traits = c(currentT,current_exT_exC),
      h2_z = pruneddisease_category$h2_z[match(c(currentT, current_exT_exC),pruneddisease_category$traits)]
    )
    phenotype_path_update<-paste0(save_path_update,'/',currentT,'_de_',unique_category_left[C],'_initialpheno.txt')

    data.table::fwrite(phenotypename_update,file = phenotype_path_update,quote = F,sep = '\t')

    h2_vector_update<-t(as.matrix(h2_vector[1,as.character(phenotypename_update$traits)]))
    h2_vector_mat_update<-h2_vector_mat[,as.character(phenotypename_update$traits)]
    Results_full_rg_update<-Results_full_rg[as.character(phenotypename_update$traits),as.character(phenotypename_update$traits)]
    Results_full_rg_array_update<-Results_full_rg_array[as.character(phenotypename_update$traits),as.character(phenotypename_update$traits),]
    Rg_mat_z_update<-Rg_mat_z[as.character(phenotypename_update$traits), as.character(phenotypename_update$traits)]

    target_num = which(phenotypename_update$traits==currentT)
    pruning_pleioh2g_corr_single_rgzscore_cutatall_more(target_num,phenotype_path_update,save_path_update,h2_vector_update,h2_vector_mat_update,Results_full_rg_update,Results_full_rg_array_update,Rg_mat_z_update,sample_rep,tolerance,seed)

  }
}
