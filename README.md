# pleioh2g
R package used in 'Pleiotropic heritability quantifies the shared genetic variance of common diseases'

## **Definition and estimation of pleiotropic heritability (h<SUP>2</SUP><SUB>pleio</SUB>)**
See Figure 1 in manuscript:
![image](https://github.com/user-attachments/assets/f7424522-f1ba-4312-9db6-e10b2e8e47a1)


The total phenotypic variance of the target disease consists of genetic variance (G) and environmental variance (E). The genetic variance of the target disease is partitioned into a disease-specific component and a pleiotropic component. The disease-specific component is not shared with the auxiliary diseases, and the pleiotropic component consists of a linear combination of the genetic values (G1, G2, …, Gn) for auxiliary diseases 1 to n.

We define **pleiotropic heritability h<SUP>2</SUP><SUB>pleio</SUB>** as the liability-scale genetic variance (estimated from SNPs) of a target disease that is shared with a specific set of auxiliary diseases. 
## **PHBC**
We estimate h<SUP>2</SUP><SUB>pleio</SUB> from GWAS summary statistics by estimating the proportion of variance explained from an estimated genetic correlation matrix and employing a Monte-Carlo bias correction procedure to account for sampling noise in genetic correlation estimates. 

### **Installation**:
```>
install.packages("devtools")
devtools::install_github("yjzhao1004/pleioh2g")
library(pleioh2g)
```
### **Data Preparation**
Our method uses multiple-trait GWAS summary statistics as input and computes genetic correlation by cross-trait LDSC (Bulik-Sullivan et al. 2015b *Nat Genet*)
* We used 1000 Genomes Project Europeans as a reference LD panel to estimate genetic correlation. 1000G European LD-scores and HapMap 3 SNPs list are available at https://data.broadinstitute.org/alkesgroup/LDSCORE.
* Summary association statistics (LDSC format .sumstat.gz) for all diseases/traits analyzed in this study are available at https://alkesgroup.broadinstitute.org/PHBC/.

### **Steps**
### Step 1: Prepare LDSC input-format data for multiple traits.
**munge_gwas_allphenotype.R** is to transform GWAS summary statistics to prepare LDSC .sumstat.gz data for all phenotypes (This function implements GenomicSEM R package; ref. Grotzinger et al. 2019)
(If you already have LDSC .sumstat.gz data, skip it and go to step 2)
 ```
 hmp3 <- system.file("extdata", "w_hm3.snplist",package = "pleioh2g") #hmp3 is the directory path to the Hapmap3 snplist file.
gwas_dir <- c(system.file("extdata/GWAS_example","401.1_gwas.txt.gz", package = "pleioh2g"),system.file("extdata/GWAS_example","250.2_gwas.txt.gz",package = "pleioh2g"),system.file("extdata/GWAS_example","296.22_gwas.txt.gz", package = "pleioh2g")) #gwas_dir is a vector of GWAS summary statistics path
 Nsamp <- c(157206,157206,157206) #Nsamp is a vector of sample size for these summary statistics - in the same order as GWAS in gwas_dir
dir.create(file.path(system.file("extdata",package = "pleioh2g"),"GWAS_sumstat"))
output_dir <- system.file("extdata/GWAS_sumstat",package = "pleioh2g") #Directory where LDSC format .sumstats.gz will be stored.
trait_names<-c("401.1","250.2","296.22") #trait_names is a vector of phenotype name for prefix name of .sumstats.gz data - in the same order as GWAS in gwas_dir
munge_gwas_allphenotype(hmp3, gwas_dir, Nsamp, trait_names,output_dir)
```
### Step 2: Compute heritability and genetic correlation point estimates from cross-trait LDSC.
**Cal_rg_h2g_alltraits.R** is to compute h<sup>2</sup> and r<sub>g</sub> point estimates using cross-trait LDSC. (This function implements ldscr R package).
```
phenotype_path<-system.file("extdata", "phenotype_name_package_test.txt",package = "pleioh2g") #phenotype_path is the directory path to the phenotype name file - a .txt file that listed the all phenotype names for all GWAS summary statistics with header named 'traits'.
gwas_dir <- c(system.file("extdata/GWAS_example","401.1_gwas.txt.gz", package = "pleioh2g"),system.file("extdata/GWAS_example","250.2_gwas.txt.gz",package = "pleioh2g"),system.file("extdata/GWAS_example","296.22_gwas.txt.gz", package = "pleioh2g")) #gwas_dir is a vector of GWAS summary statistics path
dir.create(file.path(system.file("extdata",package = "pleioh2g"),"rg_results_test"))
save_path <- system.file("extdata", "rg_results_test",package = "pleioh2g") #save_path is the directory where .rds format rg + h2g matrix will be stored.

#replace ld_path and wld_path to your own LD-scores directory path
ld_path<-system.file("extdata", "eur_w_ld_chr",package = "pleioh2g")
wld_path<-system.file("extdata", "eur_w_ld_chr",package = "pleioh2g")

#If including quantitative traits, leave sample_prev_path and population_prev_path NULL.
sample_prev_path = system.file("extdata", "samp_prev_test.txt",package = "pleioh2g") #Directory for sample prevalence, having two columns, one with header named 'diseases',the other one with header named 'pre'.
population_prev_path = system.file("extdata", "pop_prev_test.txt",package = "pleioh2g") #Directory for population prevalence, having two columns, one with header named 'diseases',the other one with header named 'pre'.

Cal_rg_h2g_alltraits(phenotype_path, gwas_munge_dir, save_path, ld_path, wld_path, sample_prev_path, population_prev_path)
```
Through this step, you can obtain h<sup>2</sup>, h<sup>2</sup> z-scores, r<sub>g</sub>, r<sub>g</sub> z-scores, genetic covariance point estimates matrix ('Results_full_h2.rds', 'Results_full_h2Z.rds', 'Results_full_rg.rds', 'Results_full_rgz.rds', 'Results_full_gcov.rds') and liability-scale h<sup>2</sup> matrix ('Results_full_h2_lia.rds') (if you set sample_prev_path and population_prev_path) in the save_path.
### Step 3: Compute heritability and genetic correlation jackknife estimates from cross-trait LDSC.
**Cal_rg_h2g_jk_alltraits.R** is to perform genomic-block jackknife and computes h<sup>2</sup> and r<sub>g</sub> jackknife estimates using cross-trait LDSC. (This function implements ldscr R package).
* We reimplement cross-trait LDSC jackknife procedure to guarantee the standard errors of all elements in genetic correlation matrix are estimated through the same genomic jackknife blocks.
* Jackknife blocks are created by partitioning 1,217,311 HapMap 3 SNPs. 
```
hmp3 <- system.file("extdata", "w_hm3.snplist",package = "pleioh2g") #hmp3 is the directory path to the Hapmap3 snplist file.
phenotype_path<-system.file("extdata", "phenotype_name_package_test.txt",package = "pleioh2g") #phenotype_path is the directory path to the phenotype name file - a .txt file that listed the all phenotype names for all GWAS summary statistics with header named 'traits'.
gwas_dir <- c(system.file("extdata/GWAS_example","401.1_gwas.txt.gz", package = "pleioh2g"),system.file("extdata/GWAS_example","250.2_gwas.txt.gz",package = "pleioh2g"),system.file("extdata/GWAS_example","296.22_gwas.txt.gz", package = "pleioh2g")) #gwas_dir is a vector of GWAS summary statistics path
dir.create(file.path(system.file("extdata",package = "pleioh2g"),"rg_results_test"))
save_path <- system.file("extdata", "rg_results_test",package = "pleioh2g") #save_path is the directory where .rds format rg + h2g matrix will be stored.

#replace ld_path and wld_path to your own LD-scores directory path
ld_path<-system.file("extdata", "eur_w_ld_chr",package = "pleioh2g")
wld_path<-system.file("extdata", "eur_w_ld_chr",package = "pleioh2g")

#If including quantitative traits, leave sample_prev_path and population_prev_path NULL.
sample_prev_path = system.file("extdata", "samp_prev_test.txt",package = "pleioh2g") #Directory for sample prevalence, having two columns, one with header named 'diseases',the other one with header named 'pre'.
population_prev_path = system.file("extdata", "pop_prev_test.txt",package = "pleioh2g") #Directory for population prevalence, having two columns, one with header named 'diseases',the other one with header named 'pre'.

#Set the number of jackknife blocks
n_block=200

Cal_rg_h2g_jk_alltraits(n_block, hmp3, phenotype_path, gwas_munge_dir, save_path, ld_path, wld_path, sample_prev_path, population_prev_path)
```
Through this step, you can obtain h<sup>2</sup> and r<sub>g</sub> jackknife estimate array ('Results_full_h2_array.rds', 'Results_full_rg_array.rds') and liability-scale h<sup>2</sup> jackknife estimate array ('Results_full_h2_lia_array.rds') (if you set sample_prev_path and population_prev_path) in the save_path.
### Step 4: Compute pleiotropic heritability 
**pruning_pleioh2g_corr_single_rgzscore_cutatall_more.R** is to in computing h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> while performing pruning and bias correction.
* We note that h<SUP>2</SUP><SUB>pleio</SUB> is a function of both the target disease and the selected set of auxiliary diseases/traits. We use the ratio of pleiotropic heritability vs. total heritability (h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP>) to quantify the proportion of genetic variance that is pleiotropic.
```
# First to determine which disease in your list is the target disease
G = 3 # Index of target disease in trait list

# Load your h2, rg, rg z-scores point estimates and jackknife estimates matrix. Take our implemented data as example: 
data("Results_full_rg_45D")
data("Results_full_rg_array_45D")
data("h2_vector_mat_45D")
data("h2_vector_45D")
data("Rg_mat_z_45D")

phenotype_path <- system.file("extdata", "Disease45_h2_Z.txt", package = "pleioh2g") # Directory path to the phenotype name and h2z score .txt file - 1st col name: 'traits'; 2nd col namd: 'h2_z'
dir.create(file.path(system.file("extdata",package = "pleioh2g"),"save_results_test"))
dir.create(file.path(system.file("extdata/save_results_test",package = "pleioh2g"),"D"))
save_path <- system.file("extdata", "save_results_test/D",package = "pleioh2g") # Directory where post-corrected pleioh2g results will be stored.
sample_rep <- 1000 # Monte Carlo sampling iterations in bias correction
tolerance <- 1e-6 # Default tolerance for binary search in bias correction
seed <- 123 # Default random seed 

pruning_pleioh2g_corr_single_rgzscore_cutatall_more(G,phenotype_path,save_path,h2_vector_45D,h2_vector_mat_45D,Results_full_rg_45D,Results_full_rg_array_45D,Rg_mat_z_45D,sample_rep,tolerance,seed)
```
You can obtain a result table with colnames and corresponding meanings:
| target_disease | target_disease_h2_est | target_disease_h2_se | selected_pleio_pheno | h2pleiotropy | h2pleiotropy_se | percentage_jackknife | percentage_jackknife_se | h2pleio_corr | h2pleio_corr_se | percentage_corr |percentage_corr_se| corrected_weight |
|----------------|-----------------------|----------------------|----------------------|--------------|-----------------|----------------------|-------------------------|--------------|-----------------|------------------|-----------------|--------------|
'target disease name' | 'h<sup>2</sup> estimate of target disease' | 'h<sup>2</sup> jackknife s.e. of target disease' | 'auxiliary phenotypes included to compute h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP>' | 'pre-correction h<SUP>2</SUP><SUB>pleio</SUB> ' | 'pre-correction h<SUP>2</SUP><SUB>pleio</SUB> jackknife s.e.' | 'pre-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<sup>2</sup>' | 'pre-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<sup>2</sup> jackknife s.e.' | 'post-correction h<SUP>2</SUP><SUB>pleio</SUB>'| 'post-correction h<SUP>2</SUP><SUB>pleio</SUB> jackknife s.e.' | 'post-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<sup>2</sup>'| 'post-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<sup>2</sup> jackknife s.e.' | 'corrected weight ξ<sub>c</sub>'|

You can also obtain a .rds data which is a vector containing each pre-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<sup>2</sup> jackknife estimate.

### Additional step: compute pleiotropic heritability in leave-category-out analyses 
**pleiotropyh2_D_T.R** is to compute h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> excluding target disease category

**pleiotropyh2_D_othersone.R** is to compute h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> excluding one other specified category

**pleiotropyh2_D_othersall.R** is to compute h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> excluding each other category one by one

**pleiotropyh2_D_T_othersone.R** is to compute h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> excluding target disease category and one other specified category

**pleiotropyh2_D_T_othersall.R** is to compute h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> excluding target disease category and each other category one by one

* We note these analyses need to use interactively dependent auxiliary diseases sets. 
    * If you perform analyses excluding target disease category or any one category, you need to first perform analyses using all auxiliary diseases and take the output result table as input ('allauxD_results_path') in pleiotropyh2_D_T.R or pleiotropyh2_D_othersone.R, because the pruning begins at the actual auxiliary diseases using in last analyses.
    * If you perform analyses excluding target disease category and one other specified category, you need to first perform analyses using all auxiliary diseases excluding target disease category and take the output result table as input ('allauxD_T_results_path') in pleiotropyh2_D_T_othersone.R, because the pruning begins at the actual auxiliary diseases using in last analyses. (See example as below)
```
#First to determine which disease in your list is the target disease
G = 3 # Index of target disease in trait list

# Load your h2, rg, rg z-scores point estimates and jackknife estimates matrix. Take our implemented data as example: 
data("Results_full_rg_45D")
data("Results_full_rg_array_45D")
data("h2_vector_mat_45D")
data("h2_vector_45D")
data("Rg_mat_z_45D")

#Before you preform leave-category-out analyses, you need to prepare a .csv or .txt phenotype file including phenotype name, h<sup>2</sup> z-scores and category information data with col names: 1st col: 'traits'; 2nd col: 'h2_z'; 3rd col: 'Category' 
#' phenotype_path <- system.file("extdata", "Disease_45d_h2Z_Category.csv",package = "pleioh2g")

dir.create(file.path(system.file("extdata",package = "pleioh2g"),"save_results_test"))
dir.create(file.path(system.file("extdata/save_results_test",package = "pleioh2g"),"D_T_others"))
save_path <- system.file("extdata", "save_results_test/D_T_others",package = "pleioh2g") # Directory where post-corrected pleioh2g results will be stored.
sample_rep <- 1000 # Monte Carlo sampling iterations in bias correction
tolerance <- 1e-6 # Default tolerance for binary search in bias correction
seed <- 123 # Default random seed 

#Specified the directory path to the result using all auxiliary disease categories excluding target disease category - the same target disease as G
allauxD_T_results_path<-system.file("extdata/save_results_test/D_T", "318_h2pleio_dresults.csv",package = "pleioh2g")

pleiotropyh2_D_T_othersall(G,phenotype_path,save_path,allauxD_T_results_path,h2_vector_45D,h2_vector_mat_45D,Results_full_rg_45D,Results_full_rg_array_45D,Rg_mat_z_45D,sample_rep,tolerance,seed)
```
You can also obtain the similar result outputs in Step 4.

