#' Prune disease selection
#'
#' @param Target_disease trait_name of target disease
#' @param trait_name trait_name of pre-prune rg_matrix
#' @param Rg_mat pre-prune rg_matrix
#' @param Rg_mat_z pre-prune rg z matrix
#' @param rg_threshold rg_threshold

#'
#' @import dplyr
#' @import data.table
#'
#' @return Rg_mat_leave
#' @export
#'
#' @examples
#'
#' selectGWAS_h2g_Z_path<-system.file("extdata", "Disease61_h2_Z.txt",
#' package = "pleioh2g")
#' selectGWAS_h2g_Z<-data.table::fread(selectGWAS_h2g_Z_path)
#' trait_name<-selectGWAS_h2g_Z$traits
#' data("Results_full_rg_61D")
#' data("Rg_mat_z_61D")
#' Target_disease<-'401.1'
#' rg_threshold<-sqrt(0.5)
#' Rg_prune<-Prune_disease_selection_DTrgzscore(Target_disease, trait_name,Results_full_rg_61D,Rg_mat_z_61D,rg_threshold)
#'
Prune_disease_selection_DTrgzscore<-function(Target_disease,trait_name,Rg_mat,Rg_mat_z,rg_threshold){
  ##select aux. exp. T
  trait_name<-setdiff(trait_name,Target_disease)
  Rg_mat_DD<-Rg_mat[as.character(trait_name),as.character(trait_name)]
  Rg_mat_DT<-Rg_mat[Target_disease,as.character(trait_name)]

  Rg_mat_z_DT<-Rg_mat_z[Target_disease,as.character(trait_name)]

  #Pruning based on genetic correction
  indices <- which((Rg_mat_DD > rg_threshold) & lower.tri(Rg_mat_DD), arr.ind = TRUE)
  result <- apply(indices, 1, function(idx) {
    paste0("(", rownames(Rg_mat_DD)[idx[1]], ", ", colnames(Rg_mat_DD)[idx[2]], ")")
  })
  result<-as.vector(result)

  if(length(result)==0){
    cat("The current rg matrix cannot be pruned at this rg_threshold.")
    Rg_mat_leave<-Rg_mat
  }else{


    pairdisease<-list()
    for(i in c(1:length(result))){
      input_string <- result[i]
      cleaned_string <- gsub("[()]", "", input_string)
      split_data <- strsplit(cleaned_string, ",\\s*")[[1]]

      if(abs(Rg_mat_z_DT[which(names(Rg_mat_z_DT)==split_data[1])])>abs(Rg_mat_z_DT[which(names(Rg_mat_z_DT)==split_data[2])])){
        pairdisease$disease1[i]<-split_data[1]
        pairdisease$rg1[i]<-Rg_mat_DT[which(names(Rg_mat_DT)==split_data[1])]
        pairdisease$disease2[i]<-split_data[2]
        pairdisease$rg2[i]<-Rg_mat_DT[which(names(Rg_mat_DT)==split_data[2])]
      }else{
        pairdisease$disease1[i]<-split_data[2]
        pairdisease$rg1[i]<-Rg_mat_DT[which(names(Rg_mat_DT)==split_data[2])]
        pairdisease$disease2[i]<-split_data[1]
        pairdisease$rg2[i]<-Rg_mat_DT[which(names(Rg_mat_DT)==split_data[1])]
      }

    }

    pairdisease<-as.data.frame(pairdisease)
    pairdisease <- pairdisease[order(-pairdisease$rg1),]

    Rg_mat_leave<-Rg_mat[setdiff(rownames(Rg_mat),unique(pairdisease$disease2)),setdiff(rownames(Rg_mat),unique(pairdisease$disease2))]

    cat("remove", paste0(unique(pairdisease$disease2),collapse = ","), "in the aux. diseases due to low rg zscore with Target disease.. \n")
  }

  return(Rg_mat_leave)

}

