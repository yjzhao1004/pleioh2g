# pleioh2g
R package used in 'Pleiotropic shared heritability quantifies the shared genetic variance of common diseases'

## **Definition and estimation of pleiotropic shared heritability (h<SUP>2</SUP><SUB>pleio</SUB>)**
See Figure 1 in manuscript:
![image](https://github.com/user-attachments/assets/f7424522-f1ba-4312-9db6-e10b2e8e47a1)


The total phenotypic variance of the target disease consists of genetic variance (G) and environmental variance (E). The genetic variance of the target disease is partitioned into a disease-specific component and a pleiotropic component. The disease-specific component is not shared with the auxiliary diseases, and the pleiotropic component consists of a linear combination of the genetic values (G<SUB>1</SUB>, G<SUB>2</SUB>, …, G<SUB>n</SUB>) for auxiliary diseases 1 to n.

We define **pleiotropic shared heritability h<SUP>2</SUP><SUB>pleio</SUB>** as the liability-scale genetic variance (estimated from SNPs) of a target disease that is shared with a specific set of auxiliary diseases. 
## **PHBC**
Our method uses multiple-trait GWAS summary statistics as input. We estimate h<SUP>2</SUP><SUB>pleio</SUB> from GWAS summary statistics by estimating the proportion of variance explained from an estimated genetic correlation matrix and employing a Monte-Carlo bias correction procedure to account for sampling noise in genetic correlation estimates. 

### **Trait Selection**
Our method recomments users to select target and auxiliary traits with heritability z-score > 6, and include auxiliary traits with r<SUB>g</SUB><SUP>2</SUP> < 0.5 with the target trait. Our software will automatically remove auxiliary traits with r<SUB>g</SUB><SUP>2</SUP> > 0.5 during computation. 



For bug reports, please email: yujiezhao@hsph.harvard.edu.


### **Installation**:
```>
install.packages("devtools")
library(devtools)
devtools::install_github("yjzhao1004/pleioh2g")
library(pleioh2g)
```
* We have published the package on CRAN (https://cran.r-project.org/web/packages/pleioh2g). You can choose to install the package from CRAN as below but without the implemented example summary statistics and 1000G EUR LD scores due to the limited memory, which cannot use "sumstats_munged_example_input()" function in Step 1. 
```>
install.packages("pleioh2g")
library(pleioh2g)
```

### **Data Preparation**
* We used 1000 Genomes Project Europeans as a reference LD panel to estimate genetic correlation. 1000G European LD-scores and HapMap 3 SNPs list are available at https://data.broadinstitute.org/alkesgroup/LDSCORE.
* Summary association statistics (LDSC format .sumstat.gz) for all diseases/traits analyzed in this study are available at https://alkesgroup.broadinstitute.org/PHBC/.

### **Steps**
### Step 1: Prepare LDSC input-format data for multiple traits.
PHBC needs LDSC .sumstats format data as input, so you first need to reformat the GWAS summary statistics as LDSC requested. We strongly recommend that you use the script munge_sumstats.py included in LDSC python package (https://github.com/bulik/ldsc) to convert your own GWAS summary statistics into LDSC format. (ref. Bulik-Sullivan et al. 2015b *Nat Genet*)
We also provide all LDSC-format .sumstat.gz data used in our analyses. (See Data preparation)
* Examples of three phenotypes ("401.1", "250.2", "296.22") have been implemented in our package. You can use codes below to reload them.

```
library(pleioh2g)
munged_sumstats = list("401.1" = sumstats_munged_example_input(example = "401.1"), "250.2" = sumstats_munged_example_input(example = "250.2"),"296.22" = sumstats_munged_example_input(example = "296.22"))
```
### Step 2: Compute pleiotropic shared heritability with bias correction
The function **pruning_pleioh2g_wrapper()** (See example as below) is to compute h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> while performing pruning and bias correction with ldsc-format GWAS summary statistics (.sumstat.gz) as input.

* We note that h<SUP>2</SUP><SUB>pleio</SUB> is a function of both the target disease and the selected set of auxiliary diseases/traits. We use the ratio of pleiotropic shared heritability vs. total heritability (h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP>) to quantify the proportion of genetic variance that is pleiotropic.
* We just use 5 jackknife blocks and 3 traits as example for quick computation test. If you set 200 jackknife blocks and include more than 50 traits, the procedure will cost more than 10 hours.

```
# Specify phenotype names
phenotype<-c("401.1","250.2","296.22")

# First to determine which disease in your list is the target disease
G = 1 # Index of target disease in trait list - this example is to compute pleiotropic shared heritability for "401.1".

# Input ldsc format .sumstat.gz data
munged_sumstats = list("401.1" = sumstats_munged_example_input(example = "401.1"), "250.2" = sumstats_munged_example_input(example = "250.2"),"296.22" = sumstats_munged_example_input(example = "296.22"))

# Specify reference LD data: ld and wld path; and hapmap 3 SNPs list
ld_path<-fs::path(fs::path_package("extdata/eur_w_ld_chr", package = "pleioh2g"))
wld_path<-fs::path(fs::path_package("extdata/eur_w_ld_chr", package = "pleioh2g"))
hmp3<-fs::path(fs::path_package("extdata/w_hm3.snplist", package = "pleioh2g"))

# If you trait is disease phenotype or the other binary trait, specify prevalence to compute the liability-scale heritability; If you don't specify this, it will compute observed-scale heritability.
sample_prev <- c(0.37,0.1,0.17)
population_prev <- c(0.37,0.1,0.17)

# Specify number of genomic-jackknife block; We use n_block = 5 as example for quick computation, but we recommand to use 200 jackknife blocks; If you don't specify this, the default number is 200.
n_block<-5
 
# Specify number of Monte Carlo sampling iterations in bias correction; If you don't specify this, the default number is 1000.
sample_rep <- 1000 

post_correction_results<-pruning_pleioh2g_wrapper(G,phenotype,munged_sumstats,ld_path, wld_path, sample_prev, population_prev,n_block, hmp3,sample_rep)
```
* If you install the package from CRAN, you should upload the 'munged_sumstats' and set the 'ld_path', 'wld_path', 'hmp3' pathway from your own file before runing the function: pruning_pleioh2g_wrapper().
```
##"gwas_munge_dir" is the directory where you store the summary statistics data. Each summary statistics file is named in the format “phenotype.sumstat.gz”, where phenotype is specified in the beginning.

target_pheno_list<-paste0(gwas_munge_dir,phenotype,'.sumstats.gz')
GWAS_list <- lapply(1:length(phenotype), function(i) fread(target_pheno_list[i]) %>% filter(!is.na(N)))
munged_sumstats = setNames(GWAS_list, phenotype)

hmp3 <- './eur_w_ld_chr/w_hm3.snplist'
ld_path<-'./eur_w_ld_chr/'
wld_path<-'./eur_w_ld_chr/'
```


An output line will provide your post-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> estimate, along with a result list `post_correction_results`, containing the following elements：
  - `target_disease` (character): The value "401.1".
  - `target_disease_h2_est` (numeric): target disease h<SUP>2</SUP>.
  - `target_disease_h2_se` (numeric): target disease h<SUP>2</SUP> s.e..
  - `selected_auxD` (character): auxiliary diseases.
  - `h2pleio_uncorr` (numeric): pre-correction h<SUP>2</SUP><SUB>pleio</SUB> estimate.
  - `h2pleio_uncorr_se` (numeric): pre-correction h<SUP>2</SUP><SUB>pleio</SUB> jackknife s.e. estimate.
  - `percentage_h2pleio_uncorr` (numeric): pre-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> estimate.
  - `percentage_h2pleio_uncorr_se` (numeric): pre-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> jackknife s.e. estimate.
  - `percentage_h2pleio_jackknife_uncorr` (numeric): vector of all pre-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> jackknife estimates in default 200 blocks.
  - `h2pleio_corr` (numeric): post-correction h<SUP>2</SUP><SUB>pleio</SUB> estimate.
  - `h2pleio_corr_se` (numeric): post-correction h<SUP>2</SUP><SUB>pleio</SUB> estimate s.e..
  - `percentage_h2pleio_corr` (numeric): post-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> estimate.
  - `percentage_h2pleio_corr_se` (numeric): post-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> jackknife s.e. estimate.
  - `percentage_h2pleio_corr_Z` (numeric): post-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> estimate Z score.
  - `corrected_weight` (numeric): corrected weight ξc in bias correction.



### Leave-category-out analysis instruction  
* If you want to compute pleiotropic shared heritability with respect to a subset of auxiliary disease, you can just specify the subset as the auxiliary disease and repeat step 1 and 2
* Estimating the contribution of a subset of auxiliary categories or diseases requires careful consideration due to an inherent pruning procedure. Our recommended approach is as follows. 
    * For analyses excluding the target disease category or a single specified category:
         * First, conduct the analyses using all auxiliary diseases.
         * Then, remove the auxiliary diseases within the specified category from `post_correction_results$selected_auxD`.
         * Update the `phenotype` and `munged_sumstats` input parameters with the remaining auxiliary diseases.
         * Rerun Step 2.
    * For analyses excluding the target disease category and an additional specified category:
         * First, conduct the analyses using all auxiliary diseases _excluding_ the target disease category.
         * Then, remove the auxiliary diseases within the specified category from `post_correction_results$selected_auxD`.
         * Update the `phenotype` and `munged_sumstats` input parameters with the remaining auxiliary diseases.
         * Rerun Step 2.
    * This procedure will ensure the consistency of pruning with the subsequent analyses. 
