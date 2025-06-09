#' Perform pruning in computing pleioh2g w.r.t non-ukb diseases
#'
#' @param G index of target disease.
#' @param phenotype_path Path to the phenotype name and h2z score file - 1st col: 'traits'; 2nd col: 'h2_z'
#' @param Nstartc Starting index of non-UKB diseases in phenotype names
#' @param save_path directory where post-corrected pleioh2g results will be stored.
#' @param h2_vector h2g vector for all traits - 1✖number_of_disease matrix - aligned as the order in phenotype file
#' @param h2_vector_mat h2g array from jackknife-block estimates for all traits - blocknum✖number_of_disease matrix - aligned as the order in phenotype file
#' @param Results_full_rg genetic correlation matrix.- number_of_disease✖number_of_disease matrix - aligned as the order in phenotype file
#' @param Results_full_rg_array genetic correlation jackknife-block array.- number_of_disease✖number_of_disease✖blocknum matrix - aligned as the order in phenotype file
#' @param Rg_mat_z rg z matrix
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
#' G = 2
#' data("Results_full_rg_45D")
#' data("Results_full_rg_array_45D")
#' data("h2_vector_mat_45D")
#' data("h2_vector_45D")
#' data("Rg_mat_z_45D")
#' phenotype_path <- system.file("extdata", "Disease45_h2_Z.txt",
#' package = "pleioh2g")
#' Nstartc<-16
#' save_path <- system.file("extdata","save_results_test_onlynonukb", package = "pleioh2g")
#' sample_rep <- 10
#' tolerance <- 1e-6
#' seed <- 123
#' pruning_pleioh2g_corr_single_onlynonukb_DTrgzscore_seriesmore(G,phenotype_path,Nstartc,save_path,h2_vector_45D,h2_vector_mat_45D,Results_full_rg_45D,
#' Results_full_rg_array_45D,Rg_mat_z_45D,sample_rep,tolerance,seed)

pruning_pleioh2g_corr_single_onlynonukb_DTrgzscore_seriesmore<-function(G,phenotype_path,Nstartc,save_path,h2_vector,h2_vector_mat,Results_full_rg,Results_full_rg_array,Rg_mat_z,sample_rep,tolerance,seed){
  phenotype_name<-data.table::fread(phenotype_path)
  current_D<-as.character(phenotype_name$traits[G])
  print(current_D)
  #step 1: check D - T rg < sqrt(0.5)

  if(!all(Results_full_rg[current_D,Nstartc:nrow(Results_full_rg)]<sqrt(0.5))){
    cat("target disease =", current_D, "needs pruning due to large rg in the aux. matrix \nPerform first pruning...\n")

    Rg_DT <- Results_full_rg[current_D, Nstartc:nrow(Results_full_rg)]
    Rg_DT_leave<-Rg_DT[which(Rg_DT<sqrt(0.5))]

    traitname1<-c(current_D,names(Rg_DT_leave))

    Rg_new<-Results_full_rg[as.character(traitname1), as.character(traitname1)]

    phenotype_path_update <- file.path(save_path, paste0(current_D, "_prune_pheno.txt"))
    phenotypename <- data.frame(
      traits = traitname1
    )
    data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

    h2_vector_update_1<-h2_vector[1,phenotypename$traits]
    h2_vector_mat_update_1<-h2_vector_mat[,phenotypename$traits]
    Results_full_rg_update_1<-Results_full_rg[phenotypename$traits,phenotypename$traits]
    Results_full_rg_array_update_1<-Results_full_rg_array[phenotypename$traits,phenotypename$traits,]
    target_num<-which(phenotypename$traits==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": first pruning on D_T rg finished \n")
    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

    rg_threshold<-sqrt(0.5)
  }else{
    traitname1<-c(current_D,phenotype_name$traits[Nstartc:nrow(Results_full_rg)])
    Rg_new<-Results_full_rg[as.character(traitname1), as.character(traitname1)]

    phenotype_path_update <- file.path(save_path, paste0(current_D, "_prune_pheno.txt"))
    phenotypename <- data.frame(
      traits = traitname1
    )
    data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

    h2_vector_update_1<-h2_vector[1,phenotypename$traits]
    h2_vector_mat_update_1<-h2_vector_mat[,phenotypename$traits]
    Results_full_rg_update_1<-Results_full_rg[phenotypename$traits,phenotypename$traits]
    Results_full_rg_array_update_1<-Results_full_rg_array[phenotypename$traits,phenotypename$traits,]
    target_num<-which(phenotypename$traits==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)
    rg_threshold<-sqrt(0.5)
    cat("target disease =", current_D, ": don't need first pruning on D_T rg... next \n")
    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")
  }


  #step 2: check pre-corr jk se

  if(precorrresults$percentage_h2pleio_se>0.5){
    cat("target disease =", current_D, "needs pruning due to >0.5 jackknife s.e. \nPerform second pruning...\n")

    traitname_2 <- using_pheno

    Rg_DDT<-Results_full_rg[as.character(traitname_2), as.character(traitname_2)]
    Rg_DDT_z<-Rg_mat_z[as.character(traitname_2), as.character(traitname_2)]

    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitname_2, Rg_DDT,Rg_DDT_z, rg_threshold = sqrt(0.4))

    phenotype_path_update <- file.path(save_path, paste0(current_D, "_prune_pheno.txt"))
    phenotypename <- data.frame(
      traits = rownames(Rg_prune)
    )
    data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

    phenotypename <- data.table::fread(phenotype_path_update)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename$traits), as.character(phenotypename$traits)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename$traits), as.character(phenotypename$traits), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename$traits)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename$traits)]
    target_num<-which(phenotypename$traits==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

    rg_threshold<-sqrt(0.4)
  }

  if(precorrresults$percentage_h2pleio_se>0.5){
    cat("target disease =", current_D, "needs pruning due to >0.5 jackknife s.e. \nPerform second pruning...\n")

    traitname_2 <- using_pheno

    Rg_DDT<-Results_full_rg[as.character(traitname_2), as.character(traitname_2)]
    Rg_DDT_z<-Rg_mat_z[as.character(traitname_2), as.character(traitname_2)]

    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitname_2, Rg_DDT,Rg_DDT_z, rg_threshold = sqrt(0.3))

    phenotype_path_update <- file.path(save_path, paste0(current_D, "_prune_pheno.txt"))
    phenotypename <- data.frame(
      traits = rownames(Rg_prune)
    )
    data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

    phenotypename <- data.table::fread(phenotype_path_update)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename$traits), as.character(phenotypename$traits)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename$traits), as.character(phenotypename$traits), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename$traits)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename$traits)]
    target_num<-which(phenotypename$traits==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

    rg_threshold<-sqrt(0.3)
  }

  if(precorrresults$percentage_h2pleio_se>0.5){
    cat("target disease =", current_D, "needs pruning due to >0.5 jackknife s.e. \nPerform second pruning...\n")

    traitname_2 <- using_pheno

    Rg_DDT<-Results_full_rg[as.character(traitname_2), as.character(traitname_2)]
    Rg_DDT_z<-Rg_mat_z[as.character(traitname_2), as.character(traitname_2)]

    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitname_2, Rg_DDT,Rg_DDT_z, rg_threshold = sqrt(0.2))

    phenotype_path_update <- file.path(save_path, paste0(current_D, "_prune_pheno.txt"))
    phenotypename <- data.frame(
      traits = rownames(Rg_prune)
    )
    data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

    phenotypename <- data.table::fread(phenotype_path_update)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename$traits), as.character(phenotypename$traits)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename$traits), as.character(phenotypename$traits), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename$traits)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename$traits)]
    target_num<-which(phenotypename$traits==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

    rg_threshold<-sqrt(0.2)
  }

  if(precorrresults$percentage_h2pleio_se>0.5){
    cat("target disease =", current_D, "needs pruning due to >0.5 jackknife s.e. \nPerform second pruning...\n")

    traitname_2 <- using_pheno

    Rg_DDT<-Results_full_rg[as.character(traitname_2), as.character(traitname_2)]
    Rg_DDT_z<-Rg_mat_z[as.character(traitname_2), as.character(traitname_2)]

    Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), traitname_2, Rg_DDT,Rg_DDT_z, rg_threshold = sqrt(0.1))

    phenotype_path_update <- file.path(save_path, paste0(current_D, "_prune_pheno.txt"))
    phenotypename <- data.frame(
      traits = rownames(Rg_prune)
    )
    data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

    phenotypename <- data.table::fread(phenotype_path_update)

    Results_full_rg_update_1 <- Results_full_rg[as.character(phenotypename$traits), as.character(phenotypename$traits)]
    Results_full_rg_array_update_1 <- Results_full_rg_array[as.character(phenotypename$traits), as.character(phenotypename$traits), ]
    h2_vector_update_1 <- h2_vector[1, as.character(phenotypename$traits)]
    h2_vector_mat_update_1 <- h2_vector_mat[, as.character(phenotypename$traits)]
    target_num<-which(phenotypename$traits==current_D)
    precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_1,
                                                        h2_vector_mat_update_1,Results_full_rg_update_1,Results_full_rg_array_update_1)
    using_pheno<-rownames(Results_full_rg_update_1)

    cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

    rg_threshold<-sqrt(0.1)
  }
  if (precorrresults$percentage_h2pleio_se>0.5){
    stop("we can't! sorry!")
  }else{


    cat("target disease =", current_D, ": perform bias correction...")

    #step 3:  check corrected weights and sampling

    phenotypename<-data.table::fread(phenotype_path_update)
    Results_full_rg_update_2 <- Results_full_rg[as.character(phenotypename$traits), as.character(phenotypename$traits)]
    Results_full_rg_array_update_2 <- Results_full_rg_array[as.character(phenotypename$traits), as.character(phenotypename$traits), ]
    h2_vector_update_2 <- h2_vector[1, as.character(phenotypename$traits)]
    h2_vector_mat_update_2 <- h2_vector_mat[, as.character(phenotypename$traits)]


    target_num<-which(phenotypename$traits==current_D)
    result <- tryCatch({
      if (rg_threshold > sqrt(0.1)) {
        pleiotropyh2_cor_computing_single(target_num, phenotype_path_update, save_path, h2_vector_update_2, h2_vector_mat_update_2,
                                          Results_full_rg_update_2, Results_full_rg_array_update_2, sample_rep, tolerance, seed)
      } else {
        pleiotropyh2_cor_computing_single_prune(target_num, phenotype_path_update, save_path, h2_vector_update_2, h2_vector_mat_update_2,
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

      result <- tryCatch({

        traitname_3 <- setdiff(phenotypename$traits, as.character(current_D))

        Rg_DDT3<-Results_full_rg[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]
        Rg_DDT_z3<-Rg_mat_z[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]

        traitnames3_all<-rownames(Rg_DDT3)

        Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), as.character(traitnames3_all), Rg_DDT3,Rg_DDT_z3, rg_threshold =new_rg_threshold)

        phenotype_path_update <- file.path(save_path, paste0(as.character(current_D), "_prune_pheno.txt"))
        phenotypename <- data.frame(
          traits = rownames(Rg_prune)
        )
        data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

        phenotypename <- data.table::fread(phenotype_path_update)

        Results_full_rg_update_3 <- Results_full_rg[as.character(phenotypename$traits), as.character(phenotypename$traits)]
        Results_full_rg_array_update_3 <- Results_full_rg_array[as.character(phenotypename$traits), as.character(phenotypename$traits), ]
        h2_vector_update_3 <- h2_vector[1, as.character(phenotypename$traits)]
        h2_vector_mat_update_3 <- h2_vector_mat[, as.character(phenotypename$traits)]

        target_num <- which(phenotypename$traits == current_D)
        rg_threshold <- new_rg_threshold
        precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_3,
                                                            h2_vector_mat_update_3,Results_full_rg_update_3,Results_full_rg_array_update_3)

        cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

        if (new_rg_threshold == sqrt(0.1)) {
          pleiotropyh2_cor_computing_single_prune(target_num, phenotype_path_update, save_path, h2_vector_update_3, h2_vector_mat_update_3,
                                                  Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
        } else {
          pleiotropyh2_cor_computing_single(target_num, phenotype_path_update, save_path, h2_vector_update_3, h2_vector_mat_update_3,
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
        result <- tryCatch({
          traitname_3 <- setdiff(phenotypename$traits, as.character(current_D))

          Rg_DDT3<-Results_full_rg[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]
          Rg_DDT_z3<-Rg_mat_z[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]

          traitnames3_all<-rownames(Rg_DDT3)

          Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), as.character(traitnames3_all), Rg_DDT3,Rg_DDT_z3, rg_threshold =new_rg_threshold)

          phenotype_path_update <- file.path(save_path, paste0(as.character(current_D), "_prune_pheno.txt"))
          phenotypename <- data.frame(
            traits = rownames(Rg_prune)
          )
          data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

          phenotypename <- data.table::fread(phenotype_path_update)

          Results_full_rg_update_3 <- Results_full_rg[as.character(phenotypename$traits), as.character(phenotypename$traits)]
          Results_full_rg_array_update_3 <- Results_full_rg_array[as.character(phenotypename$traits), as.character(phenotypename$traits), ]
          h2_vector_update_3 <- h2_vector[1, as.character(phenotypename$traits)]
          h2_vector_mat_update_3 <- h2_vector_mat[, as.character(phenotypename$traits)]

          precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_3,
                                                              h2_vector_mat_update_3,Results_full_rg_update_3,Results_full_rg_array_update_3)

          cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

          target_num <- which(phenotypename$traits == current_D)
          rg_threshold <- new_rg_threshold

          if (new_rg_threshold == sqrt(0.1)) {
            pleiotropyh2_cor_computing_single_prune(target_num, phenotype_path_update, save_path, h2_vector_update_3, h2_vector_mat_update_3,
                                                    Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
          } else {
            pleiotropyh2_cor_computing_single(target_num, phenotype_path_update, save_path, h2_vector_update_3, h2_vector_mat_update_3,
                                              Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
          }

          cat("final rg threshold: ", rg_threshold, "\n")
        },error = function(e) {
          if (rg_threshold == sqrt(0.3)) {
            new_rg_threshold <- sqrt(0.2)
          } else if (rg_threshold == sqrt(0.2)) {
            new_rg_threshold <- sqrt(0.1)
          }
          result <- tryCatch({
            traitname_3 <- setdiff(phenotypename$traits, as.character(current_D))

            Rg_DDT3<-Results_full_rg[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]
            Rg_DDT_z3<-Rg_mat_z[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]

            traitnames3_all<-rownames(Rg_DDT3)

            Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), as.character(traitnames3_all), Rg_DDT3,Rg_DDT_z3, rg_threshold =new_rg_threshold)

            phenotype_path_update <- file.path(save_path, paste0(as.character(current_D), "_prune_pheno.txt"))
            phenotypename <- data.frame(
              traits = rownames(Rg_prune)
            )
            data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

            phenotypename <- data.table::fread(phenotype_path_update)

            Results_full_rg_update_3 <- Results_full_rg[as.character(phenotypename$traits), as.character(phenotypename$traits)]
            Results_full_rg_array_update_3 <- Results_full_rg_array[as.character(phenotypename$traits), as.character(phenotypename$traits), ]
            h2_vector_update_3 <- h2_vector[1, as.character(phenotypename$traits)]
            h2_vector_mat_update_3 <- h2_vector_mat[, as.character(phenotypename$traits)]

            precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_3,
                                                                h2_vector_mat_update_3,Results_full_rg_update_3,Results_full_rg_array_update_3)

            cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

            target_num <- which(phenotypename$traits == current_D)
            rg_threshold <- new_rg_threshold

            if (new_rg_threshold == sqrt(0.1)) {
              pleiotropyh2_cor_computing_single_prune(target_num, phenotype_path_update, save_path, h2_vector_update_3, h2_vector_mat_update_3,
                                                      Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
            } else {
              pleiotropyh2_cor_computing_single(target_num, phenotype_path_update, save_path, h2_vector_update_3, h2_vector_mat_update_3,
                                                Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
            }

            cat("final rg threshold: ", rg_threshold, "\n")
          }, error = function(e) {
            cat("further+ pruning for", current_D, ":", e$message, "\n")
            traitname_3 <- setdiff(phenotypename$traits, as.character(current_D))

            Rg_DDT3<-Results_full_rg[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]
            Rg_DDT_z3<-Rg_mat_z[c(as.character(current_D),as.character(traitname_3)), c(as.character(current_D),as.character(traitname_3))]

            traitnames3_all<-rownames(Rg_DDT3)

            Rg_prune <- Prune_disease_selection_DTrgzscore(as.character(current_D), as.character(traitnames3_all), Rg_DDT3,Rg_DDT_z3, rg_threshold =sqrt(0.1))

            phenotype_path_update <- file.path(save_path, paste0(as.character(current_D), "_prune_pheno.txt"))
            phenotypename <- data.frame(
              traits = rownames(Rg_prune)
            )
            data.table::fwrite(phenotypename, phenotype_path_update, quote = FALSE, sep = "\t")

            phenotypename <- data.table::fread(phenotype_path_update)

            Results_full_rg_update_3 <- Results_full_rg[as.character(phenotypename$traits), as.character(phenotypename$traits)]
            Results_full_rg_array_update_3 <- Results_full_rg_array[as.character(phenotypename$traits), as.character(phenotypename$traits), ]
            h2_vector_update_3 <- h2_vector[1, as.character(phenotypename$traits)]
            h2_vector_mat_update_3 <- h2_vector_mat[, as.character(phenotypename$traits)]

            precorrresults<-pleiotropyh2_nocor_computing_single(target_num,phenotype_path_update,save_path,h2_vector_update_3,
                                                                h2_vector_mat_update_3,Results_full_rg_update_3,Results_full_rg_array_update_3)

            cat("target disease =", current_D, ": percentage of  pleioh2g before correction:", precorrresults$percentage_h2pleio, "; s.e. ", precorrresults$percentage_h2pleio_se,"\n")

            target_num <- which(phenotypename$traits == current_D)
            rg_threshold <- sqrt(0.1)

            pleiotropyh2_cor_computing_single_prune(target_num, phenotype_path_update, save_path, h2_vector_update_3, h2_vector_mat_update_3,
                                                    Results_full_rg_update_3, Results_full_rg_array_update_3, sample_rep, tolerance, seed)
            cat("final rg threshold: ", rg_threshold, "\n")

          })
        })
      })
    })
  }
}
