#' Perform pruning in computing pleioh2g and correct bias
#'
#' @param G index of target disease.
#' @param phenotype Vector of the phenotype name
#' @param munged_sumstats All LDSC-munged GWAS .stat.gz
#' @param ld_path Path to directory containing ld score files.
#' @param wld_path Path to directory containing weight files.
#' @param sample_prev Vector of sample prevalence, in the same order of input GWAS summary statistics.
#' @param population_prev Vector of population prevalence, in the same order of input GWAS summary statistics.
#' @param sample_rep sampling times in bias correction
#' @param tolerance tolerance for dichotomic search
#' @param seed random seed
#' @import dplyr
#' @import data.table
#' @importFrom stats sd
#' @importFrom utils write.csv
#' @return A `list` containing the following elements:
#'   - `target_disease` (character): The value "318".
#'   - `target_disease_h2_est` (numeric): target disease h2g.
#'   - `target_disease_h2_se` (numeric): target disease h2g_se.
#'   - `selected_pleio_pheno` (character): auxiliary diseases.
#'   - `h2pleio` (numeric): pleiotropic heritability estimate.
#'   - `h2pleio_se` (numeric): pleiotropic heritability jackknife s.e. estimate.
#'   - `percentage_h2pleio` (numeric): percentage of pleiotropic heritability estimate.
#'   - `percentage_h2pleio_se` (numeric): percentage of pleiotropic heritability jackknife s.e. estimate.
#'   - `percentage_h2pleio_jackknife` (numeric): vector of all percentage of pleiotropic heritability jackknife estimates.
#'   - `h2pleio_corr` (numeric): post-correction pleiotropic heritability estimate.
#'   - `h2pleio_corr_se` (numeric): post-correction pleiotropic heritability jackknife s.e. estimate.
#'   - `percentage_h2pleio_corr` (numeric): post-correction percentage of pleiotropic heritability estimate.
#'   - `percentage_h2pleio_corr_se` (numeric): post-correction percentage of pleiotropic heritability jackknife s.e. estimate.
#'   - `corrected_weight` (numeric): corrected weight in bias correction.

#' @export
#'
#' @examples

#' hmp3<-fs::path(fs::path_package("extdata/w_hm3.snplist", package = "pleioh2g"))
#' phenotype<-c("401.1","244.5","318","735.3","411.4","427.2","454.1","278.1","250.2","550.1","530.11","296.22","519.8","562.1","763")
#' munged_sumstats = list("401.1" = sumstats_munged_example_input(example = "401.1"),
#' "244.5" = sumstats_munged_example_input(example = "244.5"),
#' "318" = sumstats_munged_example_input(example = "318"),
#' "735.3" = sumstats_munged_example_input(example = "735.3"),
#' "411.4" = sumstats_munged_example_input(example = "411.4"),
#' "427.2" = sumstats_munged_example_input(example = "427.2"),
#' "454.1" = sumstats_munged_example_input(example = "454.1"),
#' "278.1" = sumstats_munged_example_input(example = "278.1"),
#' "250.2" = sumstats_munged_example_input(example = "250.2"),
#' "550.1" = sumstats_munged_example_input(example = "550.1"),
#' "530.11" = sumstats_munged_example_input(example = "530.11"),
#' "296.22" = sumstats_munged_example_input(example = "296.22"),
#' "519.8" = sumstats_munged_example_input(example = "519.8"),
#' "562.1" = sumstats_munged_example_input(example = "562.1"),
#' "763" = sumstats_munged_example_input(example = "763"))
#'
#' ld_path<-fs::path(fs::path_package("extdata/eur_w_ld_chr", package = "pleioh2g"))
#' wld_path<-fs::path(fs::path_package("extdata/eur_w_ld_chr", package = "pleioh2g"))
#' sample_prev <- c(0.37,0.1,0.1,0.1,0.1,0.1,0.1,0.14,0.1,0.1,0.17,0.14,0.24,0.1,0.10)
#' population_prev <- c(0.37,0.05,0.05,0.05,0.08,0.05,0.07,0.14,0.06,0.07,0.17,0.14,0.24,0.09,0.10)
#' n_block=3
#' G=3
#' sample_rep <- 100
#' tolerance <- 1e-6
#' seed <- 123
#' post_correction_results<-pruning_pleioh2g_wrapper(G,phenotype,munged_sumstats,ld_path, wld_path, sample_prev, population_prev,n_block, hmp3,sample_rep,tolerance,seed)


pruning_pleioh2g_wrapper<-function(G,phenotype,munged_sumstats,ld_path, wld_path, sample_prev = NULL, population_prev = NULL,n_block = 200, hmp3,sample_rep,tolerance,seed){
  phenotype_name<-as.character(phenotype)
  current_D<-phenotype_name[G]
  print(current_D)
  rg_threshold<-sqrt(0.5)
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

  #step 1: check D - T and D - D rg < sqrt(0.5)
  if(!all(Results_full_rg[lower.tri(Results_full_rg)]<rg_threshold)){
    cat("target disease =", as.character(current_D), "needs pruning due to large rg in the aux. matrix \nPerform first pruning...\n")

    Rg_DT <- Results_full_rg[as.character(current_D), ]
    if(length(setdiff(names(Rg_DT)[which(Rg_DT>rg_threshold)],as.character(current_D)))>0){
      cat("remove", paste0(setdiff(names(Rg_DT)[which(Rg_DT>rg_threshold)],as.character(current_D)),collapse = ","), "due to large rg with Target disease.. \n")
    }

    Rg_DT<-Rg_DT[which(Rg_DT<rg_threshold)]
    traitname1<-names(Rg_DT)

    Rg_DDT<-Results_full_rg[c(as.character(current_D),as.character(traitname1)), c(as.character(current_D),as.character(traitname1))]
    Rg_DDT_z<-Rg_mat_z[c(as.character(current_D),as.character(traitname1)), c(as.character(current_D),as.character(traitname1))]

    traitnames1_all<-rownames(Rg_DDT)

    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames1_all, Rg_DDT,Rg_DDT_z, rg_threshold)

    phenotypename_update<-rownames(Rg_prune)

    h2_vector_update_1<-h2_vector[1,phenotypename_update]
    h2_vector_mat_update_1<-h2_vector_mat[,phenotypename_update]
    Results_full_rg_update_1<-Results_full_rg[phenotypename_update,phenotypename_update]
    Results_full_rg_array_update_1<-Results_full_rg_array[phenotypename_update,phenotypename_update,]
    target_num<-which(phenotypename_update==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)


    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": first pruning finished \n")
    cat("target disease =", current_D, ": percentage of pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")
    rg_threshold<-sqrt(0.5)
    cat(" rg threshold: ",rg_threshold,"\n")
  }else{
    phenotypename_update<-phenotype_name
    precorrresults<-pleiotropyh2_nocor_computing_single(G,phenotypename_update,h2_vector,h2_vector_mat,Results_full_rg,Results_full_rg_array)
    using_pheno<-rownames(Results_full_rg)
    rg_threshold<-sqrt(0.5)
    cat("target disease =", current_D, ": don't need first pruning... next \n")
    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")
    cat("rg threshold: ",rg_threshold,"\n")
  }

  #step 2: check pre-corr jk se

  if(precorrresults$percentage_h2pleio_se>0.5){
    cat("target disease =", current_D, "needs pruning due to large jackknife s.e. \nPerform second pruning...\n")

    traitname_2 <- setdiff(using_pheno, as.character(current_D))

    Rg_DDT2<-Results_full_rg[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]
    Rg_DDT_z2<-Rg_mat_z[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]

    traitnames2_all<-rownames(Rg_DDT2)

    rg_threshold <- sqrt(0.4)
    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames2_all, Rg_DDT2,Rg_DDT_z2, rg_threshold = sqrt(0.4))

    phenotypename_update<-rownames(Rg_prune)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename_update)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename_update)]
    target_num<-which(phenotypename_update==current_D)

    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")
    rg_threshold<-sqrt(0.4)
    cat("rg threshold: ",rg_threshold,"\n")
  }

  if(precorrresults$percentage_h2pleio_se>0.5){
    cat("target disease =", current_D, "needs pruning due to large jackknife s.e. \nPerform third pruning...\n")

    traitname_2 <- setdiff(using_pheno, as.character(current_D))

    Rg_DDT2<-Results_full_rg[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]
    Rg_DDT_z2<-Rg_mat_z[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]

    traitnames2_all<-rownames(Rg_DDT2)

    rg_threshold <- sqrt(0.3)
    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames2_all, Rg_DDT2,Rg_DDT_z2, rg_threshold = sqrt(0.3))

    phenotypename_update<-rownames(Rg_prune)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename_update)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename_update)]
    target_num<-which(phenotypename_update==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")
    rg_threshold<-sqrt(0.3)
    cat("rg threshold: ",rg_threshold,"\n")
  }

  if(precorrresults$percentage_h2pleio_se>0.5){
    cat("target disease =", current_D, "needs pruning due to large jackknife s.e. \nPerform fourth pruning...\n")

    traitname_2 <- setdiff(using_pheno, as.character(current_D))

    Rg_DDT2<-Results_full_rg[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]
    Rg_DDT_z2<-Rg_mat_z[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]

    traitnames2_all<-rownames(Rg_DDT2)

    rg_threshold <- sqrt(0.2)
    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames2_all, Rg_DDT2,Rg_DDT_z2, rg_threshold = sqrt(0.2))

    phenotypename_update<-rownames(Rg_prune)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename_update)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename_update)]
    target_num<-which(phenotypename_update==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")
    rg_threshold<-sqrt(0.2)
    cat("rg threshold: ",rg_threshold,"\n")
  }

  if(precorrresults$percentage_h2pleio_se>0.5){
    cat("target disease =", current_D, "needs pruning due to large jackknife s.e. \nPerform fifth pruning...\n")

    traitname_2 <- setdiff(using_pheno, as.character(current_D))

    Rg_DDT2<-Results_full_rg[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]
    Rg_DDT_z2<-Rg_mat_z[c(as.character(current_D),as.character(traitname_2)), c(as.character(current_D),as.character(traitname_2))]

    traitnames2_all<-rownames(Rg_DDT2)

    rg_threshold <- sqrt(0.1)
    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitnames2_all, Rg_DDT2,Rg_DDT_z2, rg_threshold = sqrt(0.1))

    phenotypename_update<-rownames(Rg_prune)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename_update)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename_update)]
    target_num<-which(phenotypename_update==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")
    rg_threshold<-sqrt(0.1)
    cat("rg threshold: ",rg_threshold,"\n")
  }

  cat("target disease =", current_D, ": perform bias correction...")

  #step 3:  check corrected weights and sampling

  Results_full_rg_update_2 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
  Results_full_rg_array_update_2 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
  h2_vector_update_2 <- h2_vector[1, as.character(phenotypename_update)]
  h2_vector_mat_update_2 <- h2_vector_mat[, as.character(phenotypename_update)]


  target_num<-which(phenotypename_update==current_D)
  tryCatch({
    if (rg_threshold > sqrt(0.1)) {
      postcorrresults<-pleiotropyh2_cor_computing_single(target_num, phenotypename_update, h2_vector_update_2, h2_vector_mat_update_2,
                                        Results_full_rg_update_2, Results_full_rg_array_update_2, sample_rep, tolerance, seed)
    } else {
      postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotypename_update, h2_vector_update_2, h2_vector_mat_update_2,
                                              Results_full_rg_update_2, Results_full_rg_array_update_2, sample_rep, tolerance, seed)
    }
    cat("final rg threshold: ", rg_threshold, "\n")
  }, error = function(e) {
    cat("Target disease =", current_D, " needs further pruning. ", e$message, "\nPerform pruning...\n")


    if (rg_threshold == sqrt(0.5)) {
      new_rg_threshold <- sqrt(0.4)
    } else if (rg_threshold == sqrt(0.4)) {
      new_rg_threshold <- sqrt(0.3)
    } else if (rg_threshold == sqrt(0.3)) {
      new_rg_threshold <- sqrt(0.2)
    } else if (rg_threshold == sqrt(0.2)) {
      new_rg_threshold <- sqrt(0.1)
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

      cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")
      if(precorrresults$percentage_h2pleio_se>0.5){
        stop(cat("target disease =", current_D, "jackknife s.e. is too large (more than 0.5) - need to prune..."))
      }
      if (new_rg_threshold == sqrt(0.1)) {
        postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                                Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
      } else {
        postcorrresults<-pleiotropyh2_cor_computing_single(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                          Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
      }

      cat("final rg threshold: ", rg_threshold, "\n")

    }, error = function(e) {
      if (rg_threshold == sqrt(0.4)) {
        new_rg_threshold <- sqrt(0.3)
      } else if (rg_threshold == sqrt(0.3)) {
        new_rg_threshold <- sqrt(0.2)
      } else if (rg_threshold == sqrt(0.2)) {
        new_rg_threshold <- sqrt(0.1)
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

        cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

        target_num <- which(phenotypename_update == current_D)
        rg_threshold <- new_rg_threshold
        if(precorrresults$percentage_h2pleio_se>0.5){
          stop(cat("target disease =", current_D, "jackknife s.e. is too large (more than 0.5) - need to prune..."))
        }
        if (new_rg_threshold == sqrt(0.1)) {
          postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotype_path_update, save_path, h2_vector_update_3, h2_vector_mat_update_3,
                                                  Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
        } else {
          postcorrresults<-pleiotropyh2_cor_computing_single(target_num, phenotype_path_update, save_path, h2_vector_update_3, h2_vector_mat_update_3,
                                            Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
        }

        cat("final rg threshold: ", rg_threshold, "\n")
      }, error = function(e) {
        if (rg_threshold == sqrt(0.3)) {
          new_rg_threshold <- sqrt(0.2)
        } else if (rg_threshold == sqrt(0.2)) {
          new_rg_threshold <- sqrt(0.1)
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

          cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

          target_num <- which(phenotypename_update == current_D)
          rg_threshold <- new_rg_threshold
          if(precorrresults$percentage_h2pleio_se>0.5){
            stop(cat("target disease =", current_D, "jackknife s.e. is too large (more than 0.5) - need to prune..."))
          }
          if (new_rg_threshold == sqrt(0.1)) {
            postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                                    Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
          } else {
            postcorrresults<-pleiotropyh2_cor_computing_single(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                              Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
          }

          cat("final rg threshold: ", rg_threshold, "\n")
        }, error = function(e) {
          cat("further+ pruning for", current_D, ":", e$message, "\n")
          traitname_3 <- setdiff(phenotypename_update, as.character(current_D))

          Rg_DDT3<-Results_full_rg[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]
          Rg_DDT_z3<-Rg_mat_z[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]

          traitnames3_all<-rownames(Rg_DDT3)

          Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), as.character(traitnames3_all), Rg_DDT3,Rg_DDT_z3, rg_threshold = sqrt(0.1))

          phenotypename_update<-rownames(Rg_prune)

          Results_full_rg_update_3 <- Results_full_rg[as.character(phenotypename_update), as.character(phenotypename_update)]
          Results_full_rg_array_update_3 <- Results_full_rg_array[as.character(phenotypename_update), as.character(phenotypename_update), ]
          h2_vector_update_3 <- h2_vector[1, as.character(phenotypename_update)]
          h2_vector_mat_update_3 <- h2_vector_mat[, as.character(phenotypename_update)]

          precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotypename_update,h2_vector_update_3,
                                                              h2_vector_mat_update_3,Results_full_rg_update_3,Results_full_rg_array_update_3)

          cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

          target_num <- which(phenotypename_update == current_D)
          rg_threshold <- sqrt(0.1)

          postcorrresults<-pleiotropyh2_cor_computing_single_prune(target_num, phenotypename_update, h2_vector_update_3, h2_vector_mat_update_3,
                                                  Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
          cat("final rg threshold: ", rg_threshold, "\n")

        })
      })
    })
  })
  return(postcorrresults)
}
