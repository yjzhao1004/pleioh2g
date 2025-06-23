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
PHBC needs LDSC .sumstats format data as input, so you first need to reformat the GWAS summary statistics as LDSC requested. We strongly recommend that you use the script munge_sumstats.py included in LDSC python package (https://github.com/bulik/ldsc) to convert your own GWAS summary statistics into LDSC format. (ref. Bulik-Sullivan et al. 2015b *Nat Genet*)
We also provide all LDSC-format .sumstat.gz data used in our analyses. (See Data preparation)
* Examples of three phenotypes ("401.1", "250.2", "296.22") have been implemented in our package. You can use codes below to reload them.

```
library(pleioh2g)
munged_sumstats = list("401.1" = sumstats_munged_example_input(example = "401.1"), "250.2" = sumstats_munged_example_input(example = "250.2"),"296.22" = sumstats_munged_example_input(example = "296.22"))
```
### Step 2: Compute pleiotropic heritability with bias correction
**pruning_pleioh2g_wrapper.R** is to compute h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> while performing pruning and bias correction with ldsc-format GWAS summary statistics (.sumstat.gz) as input.

* We note that h<SUP>2</SUP><SUB>pleio</SUB> is a function of both the target disease and the selected set of auxiliary diseases/traits. We use the ratio of pleiotropic heritability vs. total heritability (h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP>) to quantify the proportion of genetic variance that is pleiotropic.

```
# Specify phenotype names
phenotype<-c("401.1","250.2","296.22")

# First to determine which disease in your list is the target disease
G = 1 # Index of target disease in trait list - this example is to compute pleiotropic heritability for "401.1".

# Input ldsc format .sumstat.gz data
munged_sumstats = list("401.1" = sumstats_munged_example_input(example = "401.1"), "250.2" = sumstats_munged_example_input(example = "250.2"),"296.22" = sumstats_munged_example_input(example = "296.22"))

# Specify reference LD data: ld and wld path
ld_path<-fs::path(fs::path_package("extdata/eur_w_ld_chr", package = "pleioh2g"))
wld_path<-fs::path(fs::path_package("extdata/eur_w_ld_chr", package = "pleioh2g"))

# If you trait is disease phenotype or the other binary trait, specify prevalence to compute the liability-scale heritability; If you don't specify this, it will compute observed-scale heritability.
sample_prev <- c(0.37,0.1,0.14)
population_prev <- c(0.37,0.06,0.14)

# Specify number of genomic-jackknife block; We use n_block = 5 as example for quick computation, but we recommand to use 200 jackknife blocks; If you don't specify this, the default number is 200.
n_block<-5
 
# Specify number of Monte Carlo sampling iterations in bias correction; If you don't specify this, the default number is 1000.
sample_rep <- 1000 

# Specify tolerance for binary search in bias correction; If you don't specify this, the default tolerance is 1e-6.
tolerance <- 1e-6 

# Specify random seed for Monte Carlo sampling iterations in bias correction; If you don't specify this, the default seed is 123.
seed <- 123 

post_correction_results<-pruning_pleioh2g_wrapper(G,phenotype,munged_sumstats,ld_path, wld_path, sample_prev, population_prev,n_block, hmp3,sample_rep,tolerance,seed)
```
You can obtain a result list `post_correction_results` containing the following elements:

  - `target_disease` (character): The value "401.1".
  - `target_disease_h2_est` (numeric): target disease h<SUP>2</SUP>.
  - `target_disease_h2_se` (numeric): target disease h<SUP>2</SUP> s.e..
  - `selected_auxD` (character): auxiliary diseases.
  - `h2pleio` (numeric): h<SUP>2</SUP><SUB>pleio</SUB> estimate.
  - `h2pleio_se` (numeric): h<SUP>2</SUP><SUB>pleio</SUB>jackknife s.e. estimate.
  - `percentage_h2pleio` (numeric): h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> estimate.
  - `percentage_h2pleio_se` (numeric): h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> jackknife s.e. estimate.
  - `percentage_h2pleio_jackknife` (numeric): vector of all h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> jackknife estimates in default 200 blocks.
  - `h2pleio_corr` (numeric): post-correction h<SUP>2</SUP><SUB>pleio</SUB> estimate.
  - `h2pleio_corr_se` (numeric): post-correction h<SUP>2</SUP><SUB>pleio</SUB> jackknife s.e. estimate.
  - `percentage_h2pleio_corr` (numeric): post-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> estimate.
  - `percentage_h2pleio_corr_se` (numeric): post-correction h<SUP>2</SUP><SUB>pleio</SUB> / h<SUP>2</SUP> jackknife s.e. estimate.
  - `corrected_weight` (numeric): corrected weight ξc in bias correction.



### Leave-category-out analysis instruction  

* We note that our leave-category-out analysis need to use interactively dependent auxiliary diseases sets because we want to compare the reduction of pleiotropic heritability from each category removal.     
    * If you perform analyses excluding target disease category or any one category, you need to first perform analyses using all auxiliary diseases and remove the auxiliary diseases in the specified category from `post_correction_results$selected_auxD`. Then use the remain auxiliary diseases to update the input parameters `phenotype` and `munged_sumstats` and rerun step 2, because the pruning begins at the actual auxiliary diseases using in last-step analyses.
    * If you perform analyses excluding target disease category and one other specified category, you need to first perform analyses using all auxiliary diseases excluding target disease category and remove the auxiliary diseases in the specified category from `post_correction_results$selected_auxD`.Then use the remain auxiliary diseases to update the input parameters `phenotype` and `munged_sumstats` and rerun step 2, because the pruning begins at the actual auxiliary diseases using in last-step analyses.

* If you only want to compute the pleiotropic heritability without the specified categories or pleiotropic heritability with respect to specified auxiliary diseases/traits, you can specify them in the input parameters `phenotype` and `munged_sumstats` in step 2.
