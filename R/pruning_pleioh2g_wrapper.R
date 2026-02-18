#' Perform pruning in computing pleioh2g and correct bias
#'
#' @param G index of target disease.
#' @param n_block number of jackknife blocks.
#' @param hmp3 Directory for hapmap 3 snplist.
#' @param phenotype Vector of the phenotype name
#' @param munged_sumstats All LDSC-munged GWAS .stat.gz
#' @param ld_path Path to directory containing ld score files.
#' @param wld_path Path to directory containing weight files.
#' @param sample_prev Vector of sample prevalence, in the same order of input GWAS summary statistics.
#' @param population_prev Vector of population prevalence, in the same order of input GWAS summary statistics.
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

#' @export
#'



pruning_pleioh2g_wrapper<-function(G,phenotype,munged_sumstats,ld_path, wld_path, sample_prev = NULL, population_prev = NULL,n_block = 200, hmp3,sample_rep){
  phenotype_name<-as.character(phenotype)
  current_D<-phenotype_name[G]
  message(current_D)
  rg_threshold<-0.5
  #Before computing pleioh2g, need rg and rg_jk_array
  results<-Cal_rg_h2g_alltraits(phenotype_name, munged_sumstats, ld_path, wld_path, sample_prev, population_prev)
  results_jk<-Cal_rg_h2g_jk_alltraits(n_block, hmp3, phenotype_name, munged_sumstats, ld_path, wld_path, sample_prev, population_prev)

  Results_full_rg<-results$rg
  # if input of prev is null, set default NA
  if (!is.null(sample_prev)&!is.null(population_prev)) {
    h2_vector<-results$liah2
    h2_vector_mat<-results_jk$liah2array
  }else{
    h2_vector<-results$h2
    h2_vector_mat<-results_jk$h2array
  }
  Results_full_rg_array<-results_jk$rgarray
  Rg_mat_z<-results$rgz

  #step 1: check D - T
  if(!all(Results_full_rg[as.character(current_D), -G]^2<rg_threshold)){
    message("target disease =", as.character(current_D), "Requires pruning due to large rg in the aux. matrix \nPerform first pruning...\n")

    Rg_DT <- Results_full_rg[as.character(current_D), ]
    if(length(setdiff(names(Rg_DT)[which(Rg_DT^2>rg_threshold)],as.character(current_D)))>0){
      message("remove", paste0(setdiff(names(Rg_DT)[which(Rg_DT^2>rg_threshold)],as.character(current_D)),collapse = ","), "due to large rg with target disease.. \n")
    }

    Rg_DT<-Rg_DT[which(Rg_DT^2<rg_threshold)]
    traitname1<-names(Rg_DT)

    Rg_DDT<-Results_full_rg[c(as.character(current_D),as.character(traitname1)), c(as.character(current_D),as.character(traitname1))]
    Rg_DDT_z<-Rg_mat_z[c(as.character(current_D),as.character(traitname1)), c(as.character(current_D),as.character(traitname1))]

    traitnames1_all<-rownames(Rg_DDT)


    phenotypename_update<-traitnames1_all

    h2_vector_update_1<-h2_vector[1,phenotypename_update]
    h2_vector_mat_update_1<-h2_vector_mat[,phenotypename_update]
    Results_full_rg_update_1<-Results_full_rg[phenotypename_update,phenotypename_update]
    Results_full_rg_array_update_1<-Results_full_rg_array[phenotypename_update,phenotypename_update,]
    target_num<-which(phenotypename_update==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)


    message("target disease =", current_D, ": first pruning finished \n")
    message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")
    rg_threshold<-0.5
    message(" rg threshold: ",rg_threshold,"\n")
  }else{

    phenotypename_update<-c(as.character(current_D),setdiff(phenotype_name,as.character(current_D)))
    h2_vector_update_1<-h2_vector[1,phenotypename_update]
    h2_vector_mat_update_1<-h2_vector_mat[,phenotypename_update]
    Results_full_rg_update_1<-Results_full_rg[phenotypename_update,phenotypename_update]
    Results_full_rg_array_update_1<-Results_full_rg_array[phenotypename_update,phenotypename_update,]
    target_num<-which(phenotypename_update==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)


    rg_threshold<-0.5
    message("target disease =", current_D, ": don't need first pruning... next \n")
    message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")
    message("rg threshold: ",rg_threshold,"\n")
  }

  #step 2: check pre-corr jk se
  if(precorrresults$percentage_h2pleio_uncorr_se>0.5){
    message("target disease =", current_D, "needs pruning due to large jackknife s.e. \nPerform second pruning...\n")

    traitname_2 <- setdiff(phenotypename_update, as.character(current_D))

    Rg_DDT2<-Results_full_rg[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]
    Rg_DDT_z2<-Rg_mat_z[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]

    traitnames2_all<-rownames(Rg_DDT2)

    rg_threshold <- 0.5
    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames2_all, Rg_DDT2,Rg_DDT_z2, rg_threshold = 0.5)

    phenotypename_update<-rownames(Rg_prune)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename_update)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename_update)]
    target_num<-which(phenotypename_update==current_D)

    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)

    message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")
    rg_threshold<-0.5
    message("rg threshold: ",rg_threshold,"\n")
  }
  if(precorrresults$percentage_h2pleio_uncorr_se>0.5){
    message("target disease =", current_D, "needs pruning due to large jackknife s.e. \nPerform second pruning...\n")

    traitname_2 <- setdiff(phenotypename_update, as.character(current_D))

    Rg_DDT2<-Results_full_rg[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]
    Rg_DDT_z2<-Rg_mat_z[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]

    traitnames2_all<-rownames(Rg_DDT2)

    rg_threshold <- 0.4
    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames2_all, Rg_DDT2,Rg_DDT_z2, rg_threshold = 0.4)

    phenotypename_update<-rownames(Rg_prune)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename_update)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename_update)]
    target_num<-which(phenotypename_update==current_D)

    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)

    message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")
    rg_threshold<-0.4
    message("rg threshold: ",rg_threshold,"\n")
  }

  if(precorrresults$percentage_h2pleio_uncorr_se>0.5){
    message("target disease =", current_D, "needs pruning due to large jackknife s.e. \nPerform third pruning...\n")

    traitname_2 <- setdiff(phenotypename_update, as.character(current_D))

    Rg_DDT2<-Results_full_rg[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]
    Rg_DDT_z2<-Rg_mat_z[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]

    traitnames2_all<-rownames(Rg_DDT2)

    rg_threshold <- 0.3
    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames2_all, Rg_DDT2,Rg_DDT_z2, rg_threshold = 0.3)

    phenotypename_update<-rownames(Rg_prune)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename_update)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename_update)]
    target_num<-which(phenotypename_update==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)

    message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")
    rg_threshold<-0.3
    message("rg threshold: ",rg_threshold,"\n")
  }

  if(precorrresults$percentage_h2pleio_uncorr_se>0.5){
    message("target disease =", current_D, "needs pruning due to large jackknife s.e. \nPerform fourth pruning...\n")

    traitname_2 <- setdiff(phenotypename_update, as.character(current_D))

    Rg_DDT2<-Results_full_rg[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]
    Rg_DDT_z2<-Rg_mat_z[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]

    traitnames2_all<-rownames(Rg_DDT2)

    rg_threshold <- 0.2
    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames2_all, Rg_DDT2,Rg_DDT_z2, rg_threshold = 0.2)

    phenotypename_update<-rownames(Rg_prune)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename_update)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename_update)]
    target_num<-which(phenotypename_update==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)

    message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")
    rg_threshold<-0.2
    message("rg threshold: ",rg_threshold,"\n")
  }

  if(precorrresults$percentage_h2pleio_uncorr_se>0.5){
    message("target disease =", current_D, "needs pruning due to large jackknife s.e. \nPerform fifth pruning...\n")

    traitname_2 <- setdiff(phenotypename_update, as.character(current_D))

    Rg_DDT2<-Results_full_rg[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]
    Rg_DDT_z2<-Rg_mat_z[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]

    traitnames2_all<-rownames(Rg_DDT2)

    rg_threshold <- 0.1
    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames2_all, Rg_DDT2,Rg_DDT_z2, rg_threshold = 0.1)

    phenotypename_update<-rownames(Rg_prune)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename_update)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename_update)]
    target_num<-which(phenotypename_update==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)

    message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")
    rg_threshold<-0.1
    message("rg threshold: ",rg_threshold,"\n")
  }

  message("target disease =", current_D, ": perform bias correction...")

  #step 3:  check corrected weights and sampling

  Results_full_rg_update_2 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
  Results_full_rg_array_update_2 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
  h2_vector_update_2 <- h2_vector[1, as.character(phenotypename_update)]
  h2_vector_mat_update_2 <- h2_vector_mat[, as.character(phenotypename_update)]
  target_num<-which(phenotypename_update==current_D)
  message("target trait index:",target_num,"\n")

  target_num<-which(phenotypename_update==current_D)
  tryCatch({
    if (rg_threshold > 0.1) {
      postcorrresults<-pleiotropyh2_cor_computing_single(target_num, phenotypename_update, h2_vector_update_2, h2_vector_mat_update_2,
                                        Results_full_rg_update_2, Results_full_rg_array_update_2, sample_rep)
      message("target disease =", current_D, ": ","post-correction h2pleio/h2 is ", postcorrresults$percentage_h2pleio_corr," ;s.e. ",postcorrresults$percentage_h2pleio_corr_se )

    } else {
      postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotypename_update, h2_vector_update_2, h2_vector_mat_update_2,
                                              Results_full_rg_update_2, Results_full_rg_array_update_2, sample_rep)
      message("target disease =", current_D, ": ","post-correction h2pleio/h2 is ", postcorrresults$percentage_h2pleio_corr," ;s.e. ",postcorrresults$percentage_h2pleio_corr_se )

    }
    message("Final rg threshold: ", rg_threshold, "\n")
  }, error = function(e) {
    message("Target disease =", current_D, " needs further pruning. ", e$message, "\nPerform pruning...\n")


    if (rg_threshold == 0.5) {
      new_rg_threshold <- 0.4
    } else if (rg_threshold == 0.4) {
      new_rg_threshold <- 0.3
    } else if (rg_threshold == 0.3) {
      new_rg_threshold <- 0.2
    } else if (rg_threshold == 0.2) {
      new_rg_threshold <- 0.1
    }

    tryCatch({

      traitname_3 <- setdiff(phenotypename_update, as.character(current_D))

      Rg_DDT3<-Results_full_rg[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]
      Rg_DDT_z3<-Rg_mat_z[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]

      traitnames3_all<-rownames(Rg_DDT3)

      Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), as.character(traitnames3_all), Rg_DDT3,Rg_DDT_z3, rg_threshold =new_rg_threshold)

      phenotypename_update<-rownames(Rg_prune)

      Results_full_rg_update_3 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
      Results_full_rg_array_update_3 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
      h2_vector_update_3 <- h2_vector[1, as.character(phenotypename_update)]
      h2_vector_mat_update_3 <- h2_vector_mat[, as.character(phenotypename_update)]

      target_num <- which(phenotypename_update == current_D)
      rg_threshold <- new_rg_threshold
      precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_3,
                                                          h2_vector_mat_update_3,Results_full_rg_update_3,Results_full_rg_array_update_3)
      message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")


      if(precorrresults$percentage_h2pleio_uncorr_se>0.5){
        stop(message("target disease =", current_D, "jackknife s.e. is too large (more than 0.5) - need to prune..."))
      }
      if (new_rg_threshold == 0.1) {
        postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                                Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep)
        message("target disease =", current_D, ": ","post-correction h2pleio/h2 is ", postcorrresults$percentage_h2pleio_corr," ;s.e. ",postcorrresults$percentage_h2pleio_corr_se )

      } else {
        postcorrresults<-pleiotropyh2_cor_computing_single(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                          Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep)
        message("target disease =", current_D, ": ","post-correction h2pleio/h2 is ", postcorrresults$percentage_h2pleio_corr," ;s.e. ",postcorrresults$percentage_h2pleio_corr_se )

      }

      message("final rg threshold: ", rg_threshold, "\n")

    }, error = function(e) {
      if (rg_threshold == 0.4) {
        new_rg_threshold <- 0.3
      } else if (rg_threshold == 0.3) {
        new_rg_threshold <- 0.2
      } else if (rg_threshold == 0.2) {
        new_rg_threshold <- 0.1
      }
      tryCatch({
        traitname_3 <- setdiff(phenotypename_update, as.character(current_D))

        Rg_DDT3<-Results_full_rg[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]
        Rg_DDT_z3<-Rg_mat_z[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]

        traitnames3_all<-rownames(Rg_DDT3)

        Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), as.character(traitnames3_all), Rg_DDT3,Rg_DDT_z3, rg_threshold =new_rg_threshold)

        phenotypename_update<-rownames(Rg_prune)

        Results_full_rg_update_3 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
        Results_full_rg_array_update_3 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
        h2_vector_update_3 <- h2_vector[1, as.character(phenotypename_update)]
        h2_vector_mat_update_3 <- h2_vector_mat[, as.character(phenotypename_update)]

        precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_3,
                                                            h2_vector_mat_update_3,Results_full_rg_update_3,Results_full_rg_array_update_3)

        message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")

        target_num <- which(phenotypename_update == current_D)
        rg_threshold <- new_rg_threshold
        if(precorrresults$percentage_h2pleio_uncorr_se>0.5){
          stop(message("target disease =", current_D, "jackknife s.e. is too large (more than 0.5) - need to prune..."))
        }
        if (new_rg_threshold == 0.1) {
          postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotype_path_update, h2_vector_update_3, h2_vector_mat_update_3,
                                                  Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep)
          message("target disease =", current_D, ": ","post-correction h2pleio/h2 is ", postcorrresults$percentage_h2pleio_corr," ;s.e. ",postcorrresults$percentage_h2pleio_corr_se )

        } else {
          postcorrresults<-pleiotropyh2_cor_computing_single(target_num, phenotype_path_update, h2_vector_update_3, h2_vector_mat_update_3,
                                            Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep)
          message("target disease =", current_D, ": ","post-correction h2pleio/h2 is ", postcorrresults$percentage_h2pleio_corr," ;s.e. ",postcorrresults$percentage_h2pleio_corr_se )

        }

        message("Final rg threshold: ", rg_threshold, "\n")
      }, error = function(e) {
        if (rg_threshold == 0.3) {
          new_rg_threshold <- 0.2
        } else if (rg_threshold == 0.2) {
          new_rg_threshold <- 0.1
        }
       tryCatch({
          traitname_3 <- setdiff(phenotypename_update, as.character(current_D))

          Rg_DDT3<-Results_full_rg[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]
          Rg_DDT_z3<-Rg_mat_z[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]

          traitnames3_all<-rownames(Rg_DDT3)

          Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), as.character(traitnames3_all), Rg_DDT3,Rg_DDT_z3, rg_threshold =new_rg_threshold)

          phenotypename_update<-rownames(Rg_prune)

          Results_full_rg_update_3 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
          Results_full_rg_array_update_3 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
          h2_vector_update_3 <- h2_vector[1, as.character(phenotypename_update)]
          h2_vector_mat_update_3 <- h2_vector_mat[, as.character(phenotypename_update)]

          precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_3,
                                                              h2_vector_mat_update_3,Results_full_rg_update_3,Results_full_rg_array_update_3)

          message("target disease =", current_D, ": pre-correction h2pleio/h2: ", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")

          target_num <- which(phenotypename_update == current_D)
          rg_threshold <- new_rg_threshold
          if(precorrresults$percentage_h2pleio_uncorr_se>0.5){
            stop(message("target disease =", current_D, "jackknife s.e. is too large (more than 0.5) - need to prune..."))
          }
          if (new_rg_threshold == 0.1) {
            postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                                    Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep)
            message("target disease =", current_D, ": ","post-correction h2pleio/h2 is ", postcorrresults$percentage_h2pleio_corr," ;s.e. ",postcorrresults$percentage_h2pleio_corr_se )

          } else {
            postcorrresults<-pleiotropyh2_cor_computing_single(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                              Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep)
            message("target disease =", current_D, ": ","post-correction h2pleio/h2 is ", postcorrresults$percentage_h2pleio_corr," ;s.e. ",postcorrresults$percentage_h2pleio_corr_se )

          }

          message("Final rg threshold: ", rg_threshold, "\n")
        }, error = function(e) {
          message("Further+ pruning for", current_D, ":", e$message, "\n")
          traitname_3 <- setdiff(phenotypename_update, as.character(current_D))

          Rg_DDT3<-Results_full_rg[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]
          Rg_DDT_z3<-Rg_mat_z[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]

          traitnames3_all<-rownames(Rg_DDT3)

          Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), as.character(traitnames3_all), Rg_DDT3,Rg_DDT_z3, rg_threshold = 0.1)

          phenotypename_update<-rownames(Rg_prune)

          Results_full_rg_update_3 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
          Results_full_rg_array_update_3 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
          h2_vector_update_3 <- h2_vector[1, as.character(phenotypename_update)]
          h2_vector_mat_update_3 <- h2_vector_mat[, as.character(phenotypename_update)]

          precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_3,
                                                              h2_vector_mat_update_3,Results_full_rg_update_3,Results_full_rg_array_update_3)

          message("target disease =", current_D, ": pre-correction h2pleio/h2:", precorrresults$percentage_h2pleio_uncorr, "; s.e. ", precorrresults$percentage_h2pleio_uncorr_se,"\n")

          target_num <- which(phenotypename_update == current_D)
          rg_threshold <- 0.1

          postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                                  Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep)
          message("Final rg threshold: ", rg_threshold, "\n")
          message("target disease =", current_D, ": ","post-correction h2pleio/h2 is ", postcorrresults$percentage_h2pleio_corr," ;s.e. ",postcorrresults$percentage_h2pleio_corr_se )

        })
      })
    })
  })

  return(postcorrresults)
}
