# Beta diversity change null modelling workflow for high perfomance computor clusters

This repository contains R and Slurm workflow for calculating observed and 
null model standardized alpha diversity, beta diversity, and local contribution 
to beta diversity (LCBD). Standardized effect sizes (SES) are calculated for two species pools 
(contemporary and native only), that represent current and historical communities,
as well as for the change in diversity between species pools.

This workflow is designed to run using R on local machines with 
the more expensive calculations ran using high performance computing clusters 
via the slurm interface and shell scripts.

The workflow was developed for two analyses

The goal of this repository is 1) provide code for these analyses, and
2) serve as general workflow


For detailed methodology and justification see:

>[BLANK](BLANK)

>[BLANK](BLANK)



If you use or adapt this workflow, please cite:

[Manuscript citation]

[Manuscript citation]

[Workflow citation / DOI]


For questions about this analysis, please contact:

```bash
Name: William K. Annis

Email: wannis@fsu.edu, williamkannis@gmail.com

OrcID: 0009-0003-3541-8503
```

This repository serves as a guideline for performing beta diversity change
null model analyses on HPC clusters. To reproduce the results of the associated
manuscripts, see the following repositories:

[Drivers of multifaceted beta diversity change in invaded stream fish communities](BLANK)

[Syndromes of multidimensional beta diversity change in invaded metacommunities](BLANK)

## Workflow overview

```text
                            Input data
                                │
                                ▼
                        00_prepare_inputs.R
                                │
                                ▼
                        Upload to HPC storage
                                │
                   ┌───────────────────────────┐
                   │                           │
                   ▼                           ▼
        ┌──────────────────────────────────────────────────┐
        │                  HPC / Slurm                     │
        │                                                  │
        │  Observed diversity       Null model iterations  │
        │          │                           │           │
        │          └───────────────────────────┘           │ 
        │                       │                          │
        │                       ▼                          │
        │                 Prepare SES inputs               │
        │                       │                          │
        │                       ▼                          │
        │               Calculate effect sizes             │
        └───────────────────────┬──────────────────────────┘
                                │
                                ▼
                        Download HPC results
                                │
                                ▼
                          15_ses_comp.R
                                │
                                ▼
                        Final effect-size data
```



## Workflow end products


| Diversity facet | Alpha | Beta | LCBD | Change |
| --------------- | ----: | ---: | ---: | -----: |
| Taxonomic       |     — |    ✓ |    ✓ |      ✓ |
| Functional      |     ✓ |    ✓ |    ✓ |      ✓ |
| Phylogenetic    |     ✓ |    ✓ |    ✓ |      ✓ |




## Applying workflow to new data
MENTION THAT ONLY STEP 1 IS SPECIFIC TO OUR DATA. FUTURE RELEASES WILL HAVE MORE
GENERAL NULL MODEL ALGORITHMS

MENTION THAT THIS FOCUSED ON PUEDO HISTORICAL APPROACHES WITH TWO SPECIES POOL
MAKE MORE GENERAL FUNCTIONS TO HANDLE OTHER TYPES

## Required Software
**R version**: 4.5.0

R packages

* ```'ade4'```version: 1.7.23
* ```'adespatial'``` version: 0.3.28
* ```'ape'``` version: 5.8.1
* ```'BAT'``` version: 2.11.0
* ```'DescTools'``` version: 0.99.60
* ```'dplyr'``` version: 1.1.4
* ```'ggplot2'``` version: 4.0.2
* ```'matrixStats'``` version: 1.5.0
* ```'parallel'``` version: 4.5.0
* ```'purrr'``` version: 1.1.0
* ```'tibble'``` version: 3.2.1
* ```'VGAM'```  version: 1.1.14



## Input data
First, download the entire repository. This will result in the following
file structure for uses to place  diversity input data (e.g. community, trait, 
phylogeny), formatted high performance computation input and intermediate data, 
and the resulting diversity outputs.

```bash
│
├── Diversity Input Data
│   ├── mod_com_diversity_input.rds
│   ├── his_com_diversity_input.rds
│   ├── trait_diversity_input.rds
│   └── phylo_tree.rds
│
├── HPC
│   │── scripts
│   │── beta_obs_input_data
│   ├── beta_null_input_data
│   ├── alpha_null_input_data
│   │── obs_out
│   │── null_out
│   │── ses_inputs
│   │── ses_outputs
│   └── scripts_do_not_run
│
└── Diversity Output Data

```
Place community data for both species pools (contemporary - "mod" and native - 
"his"), trait data, and phylongetic trees in ```Diversity Input Data``` with the 
following names: 

* ```mod_com_diversity_input.rds```: Community data for contemporary species 
pool. Data.frame or matrix with rows for sites and columns for species. Can 
contain an optional column ```HUC_12``` which represent regions and can be used 
to define regional species pools.
* ```his_com_diversity_input.rds```: Community data for native-only species pool. 
Same structure as ```mod_com_diversity_input.rds```
* ```trait_diversity_input.rds```: Trait data for all species in community data. 
Data.frame or matrix with rows for species and columns for traits.
* ```phylo_tree.rds```: Phylogenetic tree for all species in community data

The remaining directories will populate with data throughout the workflow.

## Workflow

Scripts needing to be ran by the use can be found in ```hpc/scripts``` 
directory and are number based on order of workflow. Scripts found
in ```hpc/scripts_do_not_run``` are either scripts ran on the HPC cluster using 
shell scripts, or contain functions used by other scripts. **DO NOT** directly run 
scripts in ```hpc/scripts_do_not_run```.

Completing the full workflow will results in observed values and null model 
standardized effect sizes for taxonomic, functional, and phylogenetic beta 
diversity/LCBD at two time steps, and change between those time steps. 
Additionally, functional and phylogenetic alpha diversity at the first time step
will be calculated in observed and effect size. Users can choose what metrics
that want to calculate by running only the scripts that correspond to the metrics
of interest. See table below for script naming codes:

```bash
tax = taxonomic
fun = functional
phy = phylogenetic
```
Additionally, users can modify script to incorporate more or less time periods,
of other diversity metrics of interest.

### 1. Prepare input data - perform on local machine
**Script:** ```00_null_input_creation.R```

<ins>Purpose:</ins> Prepares input data for the calculation of observed and null
diversity values on high performance computer clusters. By default, 999 iteration
of random community matrices (regionally constrained and taxonomic beta null 
models), random functional-spaces/phylogenetic-trees (taxa-swap null models).
Users can change the number of null iterations created by modifying the 
following code: ```max_iter  <- 999```.

Input data are prepared in chunks of random iterations to run on multiple PC nodes
simultaneously, with each chunk being ran on one node. Each node will be able
to run multiple null model iterations simultaneously across a user specified
amount of cores. By separating iterations across nodes, users can use more 
computer cores simultaneously, reduce memory pressure, reduce HPC queue times,
and more efficiency use HPC resources.

We offer recommended chunk sizes in the script based on diversity calculations
used in the studies [here](BLANK) and [here](BLANK).
For example, the 1998 iterations (999 iterations each pool) of functional beta 
diversity we divided into 999 chunks with 2 iteration per node due to the high
memory requirements of kernel density based functional metrics. COnversely,
the 1998 iteration of phylogenetic beta diversity were divided into 6 chunks
with 333 iterations per chunk

Users should determine chunk sizes based their own computational needs, and
can choose the number of chunks (i.e., nodes) per species pool by altering the 
following code:

```bash
# How many cores per nodes to get 999 total
sapply(c(1, 3, 9, 27, 37, 111),function(x) 999/x)

# Create splits based on number of nodes needed
n_groups_t <-      # taxonomic
n_groups_f <- 500  # functions
n_groups_p <- 3    # phylogenetic
```

<ins>NOTE:</ins> This script only contains code to shuffle communities using 
taxa-swap and a regionally-constrained taxa-swap null model algorithms, which 
are ran using functions called from ```null_model_algorithms.R```. While the 
algorithms provided were best suited for functional and phylogenetic beta 
diversity, they may not be suited for all null model purposes. Users should 
research the best model for their usage and modify ```00_null_input_creation.R``` 
and ```null_model_algorithms.R``` accordingly.

<ins>Outputs:</ins> Lists of input data for each metric and processing chunk.
For all metrics but taxonomic, list contains community and 
functional/phylogenetic data. Depending on null model, either the community or
functional/phylogenetic data are a list of 999 random iterations.


### 2. Upload entire ```HPC``` directory to high performance cluster storage.
For all steps utilizing HPC clusters, users will run shell scripts that will
run their respective R script using specified HPC resources. Below we only
list the shell scripts, but R scripts can be found in ```scripts_do_not_run``` and
have the same name as their shell scripts. 

We provide shell scripts with recommended Slurm arguments, but these may need to 
be altered based on the computational requirements of the user's data.

See below for an example of commonly used Slurm arguments:

```bash
# Process using 6 HPC nodes
#SBATCH --array=1-6  

# Use 37 cpu cores per node
#SBATCH --cpus-per-task=37

# Use 45gb of ram per node
#SBATCH --mem=45gb

# 6 hours of walltime
#SBATCH --time=6:00:00
```



### 3 .Calculate observed diversity values - perform using HPC
**Scripts:** ```01_obs_tax_beta.sh```, ```02_obs_fun_beta.sh```
, ```03_obs_phy_beta.sh```, ```04_obs_fun_alpha.sh``` , ```05_obs_phy_alpha.sh```

<ins>Purpose:</ins> Calculates the observed beta diversity for the contemporary 
and native only species pools, and the observed alpha diversity of the native 
only species pool. Uses one high performance computer nodes for each time step. 

<ins>Outputs:</ins>

* [Observed beta diversity data](#observed-beta-diversity-data)
* [Observed alpha diversity data](#observed-alpha-diversity-data)

### 4. Calculate null diversity values  - perform using HPC
**Scripts:** ```06_null_tax_beta.sh```, ```07_null_fun_beta.sh```
, ```08_null_phy_beta.sh```, ```09_null_fun_alpha.sh```, ```10_null_phy_alpha.sh```

<ins>Purpose:</ins> Calculates the null iterations for beta diversity in the 
contemporary and native only species pools, and null iterations for alpha 
diversity in the native only species pool. 
Set number of nodes based on number of chunks for each diversity metric. For
example, if the functional beta null input data is separated into 1000 chunks,
then ```#SBATCH --array=1-1000```. Set memory to minimum value that allows
job to run. For example: ```#SBATCH --mem=18gb```

<ins>DO NOT</ins> set cpu cores to value higher than number of iterations in 
each chunk, as this will waste resources. For example, if 2 iterations per chunk,
the maximum allowable number of cores is ```SBATCH --cpus-per-task=2 ```.

<ins>TIP:</ins> We recommend that users experiment with memory and CPU 
requirements with smaller number of null iterations before running full job.

<ins>Outputs:</ins> [Null iterations](#null-iterations)


### 5. Prepare effect size input data - perform using HPC
**Scripts:** ```11_tax_beta_null_model_prep.sh```, ```12_beta_null_model_prep.sh```
, ```13_alpha_null_model_prep.sh```

<ins>Purpose:</ins> Calculates difference in LCBD between contemporary and 
native pools (delta) for the observed values, and for each null iteration. 
Consolidates outputs into single files, with separate delta, native, and 
contemporary species pool values.

<ins>Outputs:</ins> [Intermediate files used for effect size calculations](#formatted-ses-inputs)


### 6. Calculate effect sizes - perform using HPC
**Script:**```14_batch_ses.sh```

<ins>Purpose:</ins> Estimates null model effect sizes of each metric, as well
as null model diagnostics. See [here](#effect-size-calculations) 
for more information on effect size calculations.

<ins>TIP:</ins> Adjust ```#SBATCH --array=1-26``` shell 
script argument to reflect the number of diversity metrics of interest 
(i.e., number of items in ```ses_inputs``` directory).

<ins>Outputs:</ins> [Summarized null model outputs](#summarized-null-model-outputs)

### 7. Download the entire ```HPC``` directory to local machine.


### 8. Summarize null model results - perform on local machine
**Scripts:** ```15_ses_comp.R```

<ins>Purpose:</ins> Compiles and formats the resulting SES, ES, and diagnostic 
stats across files. Also visualizes normality diagnostics to allow users to 
decide between SES or ES values for further analyses.

<ins>Outputs:</ins> 

* [```delta_lcbd.rds```](#effect-size-dataframes), 
* [```native_lcbd.rds```](#effect-size-dataframes), 
* [```native_alpha.rds```](#effect-size-dataframes)


### Helper function scripts 
```scripts_do_not_run```

* ```null_model_algorithms.R```: Contains algorithms to randomize community, 
trait, or phylogenetic data for null model analysis
* ```diversity_batch_functions.R```: Contains functions that estimate multiple 
iterations of diversity metrics using parallel computation.
* ```null_model_effect_size_function.R```: Function to summarize null 
distributions and calculates standardize effect sizes in the traditional z-score 
method (SES), empirical p-values, and p-value based effect sizes (ES), and 
reports optional diagnostic metrics.

## Null model methods

### Taxonomic beta diveristy change

### Phylogenetic and functional beta diversity

### Phylogenetic and functional alpha diversity


## Effect size calculations
```null_model_effect_size_function.R```

This workflows estimates null model effect sizes in two differnt ways:

1. Z-score based standardized effect sizes (SES)

$$
\text{SES} = \frac{\text{obs}_{\text{mean}} - \text{null}_{\text{mean}}}{\text{null}_{\text{sd}}}
$$

2. Empirical p-value based effect sizes (ES)

$$
\text{ES} = \text{probit}(1 - p) \quad \text{where} \quad p = \frac{r + 1}{n + 1}
$$

> **1. Z-score based standardized effect sizes (SES)**
> 
> $\text{SES} = \frac{\text{obs}_{\text{mean}} - \text{null}_{\text{mean}}}{\text{null}_{\text{sd}}}$

> **2. Empirical p-value based effect sizes (ES)**
> 
> $\text{ES} = \text{probit}(1 - p) \quad \text{where} \quad p = \frac{r + 1}{n + 1}$

$$
\text{1. Z-score based standardized effect sizes (SES)}
$$
$$
\text{SES} = \frac{\text{obs}_{\text{mean}} - \text{null}_{\text{mean}}}{\text{null}_{\text{sd}}}
$$

$$
\text{2. Empirical p-value based effect sizes (ES)}
$$
$$
\text{ES} = \text{probit}(1 - p) \quad \text{where} \quad p = \frac{r + 1}{n + 1}
$$


This function summarizes null distributions and calculates 
standardize effect sizes in the traditional z score method (SES), empirical 
p-values, and p-value based effect sizes (ES), and reports optional diagnostic 
metrics used to select between the two effect size methods. Asymmetrical null 
distributions should be assessed using empirical p-value based effect sizes 
rather than z-score based SES. See 
[Botta-Dukát (2018)](https://doi.org/10.1556/168.2018.19.1.8) for more 
information on selecting SES or p-value based ES.
## Output data strucutre
See below for detailed descriptions of intermidate and final data products from
this workflow.

### Observed beta diversity data 
```HPC/obs_out```

**Files:**```his_fun_beta_obs.rds```, ```his_phy_beta_obs.rds```
, ```his_tax_beta_obs.rds```, ```mod_fun_beta_obs.rds```, ```mod_phy_beta_obs.rds```
, ```mod_tax_beta_obs.rds```

Files contain lists of pairwise beta diversity distance objects and local 
contribution to beta diversity (LCBD) data frames. Files contain values for the 
contemporary (mod) and native (his) species pools, and taxonomic (tax) 
functional (fun), and phylogenetic (phy) diversity facets. All files have the 
following structure:

* ```$beta```: list of 3 containing Sorensen pairwise beta diversity distance 
objects
    * ```$Btotal```: distance object of total
    beta diversity
    * ```$Brepl```: replacement beta diversity component
    *  ```$Brich```: richness difference beta diversity component
*  ```$LCBD```: Data frame of 3 columns containing LCBD values. Column names are 
specific to diveristy facets, and these differences are denoted with *X*:
    * ```X_Btotal```: LCBD of total beta diversity in facet *X*
    * ```X_Brepl```: LCBD of replacement component in facet *X*
    * ```X_Brich```: LCBD of richness difference component in facet *X*

### Observed alpha diversity data 
```HPC/obs_out```

**Files:**```his_fun_alpha_obs.rds```, ```his_phy_alpha_obs```

Files contain named numeric objects with observed functional richness 
(volume of kernel density hypervolume) or phylogenetic richness (number of 
branches in phylogenetic tree). 

### Null iterations 
```HPC/null_out```

**Files:** Each chunk has its own file name in the following format:

> *pool*\_*facet*\_*metric*\_null\_*iteration*.rds

* *pool*: species pool (contemporary = mod; native = his)
* *facet*: diversity facet (taxonomic = tax; functional = fun; phylogenetic = phy)
* *metric*: alpha or beta diversity
* *iteration*: range of iterations included in file. e.g., 001-002

These data are the result of the randomization of traits and phylogenies (beta) 
or randomization of communities (alpha). These data will be compiled and used to 
create null distributions used to create null model standardized diversity values.

Null model iterations are divided into separate files for alpha and beta 
diversity based on diversity facet. Within diversity facets, iterations are 
broken into chunks based on the number of HPC nodes were used to create the 
data.

All files contain a list with items for each null iteration. Each beta diversity 
iteration contains a list with the same structure as 
[observed beta diversity](#observed-beta-diversity-data) and each alpha 
diversity iteration contains a named numeric object with the same structure as 
[observed alpha diversity](#observed-alpha-diversity-data).

### Formatted ses inputs 
```HPC/ses_input```

**Files:** Each diversity metricand species pool has its own file name in the following format:

> *facet*\_*pool*\_*metric*\_ses_input.rds

* *facet*: diversity facet (taxonomic = tax; functional = fun; phylogenetic = phy)
* *pool*: species pool (contemporary = mod; native = his; change = d)
* *metric*: diversity metric (Btotal = total beta diversity; Brepl = replacement 
component; Bric = richness difference component; local contribution to beta 
diversity = LCBD; alpha = alpha diversity)

These data are intermediate data product used to streamline the calculation of
standardized effect sizes.

All files contains a list with the following structure:

* ```$obs```: observed diversity values of that metric and pool
* ```$null_list: list of null iterations
   
### Summarized null model outputs 
```HPC/ses_output```

**Files:** Each diversity metric and species pool has its own file name in the following format:

> *facet*\_*pool*\_*metric*\_ses_out.rds

* *facet*: diversity facet (taxonomic = tax; functional = fun; phylogenetic = phy)
* *pool*: species pool (contemporary = mod; native = his; change = d)
* *metric*: diversity metric (Btotal = total beta diversity; Brepl = replacement 
component; Bric = richness difference component; local contribution to beta 
diversity = LCBD; alpha = alpha diversity)

These data are used to evaluate the properties of the null distributions to 
choose between standardized effect sizes and empirical effect sizes, and can be 
merged together for plotting and to create the final data.frames for the use in analyses.

All files contain a list with the following structure:

* ```$obs```: observed diversity values
* ```$null_mean```: means of null distributions
* ```$null_sd```: standards deviation of null distributions
* ```$ses```: z-score based standardized effect sizes: ```(obs - null_mean)/null_sd)```
* ```$empirical_pvalue```: proportion of null distribution with higher values than observed value
* ```$empirical_es```: Empirical p-value based effect sizes ```probit(1-empirical_pvalue)```
* ```$skew```: skewness of null distribution
* ```$kurt```: kurtosis of null distribution

The structure of the values in each list item are dependent on diversity metric:

* beta diversity components (Btotal, Brepl, and Brich): distance objects,
* LCBD: data frame with rows for each site and a column for LCBD of Btotal, Brepl, and Bric 
* alpha: named numeric objects.


### Effect size dataframes 
```Diversity Output Data```

**Files:** ```delta_lcbd.rds```, ```native_lcbd.rds```, ```native_alpha.rds```

The combined observed and summarized effect size data for native alpha, native 
LCBD, and delta LCBD. This is just one example of the type of summarized data
that can result from this workflow. These data could be used in hypothetical
analyses of beta diversity change.


* ```delta_lcbd.rds```: Observed and null model empirical p-value based effect 
sizes (ES) of taxonomic (no ES), functional and phylogenetic changes in local 
contributions to beta diversity (LCBD) between the contemporary and native species pools. 
 Data frame consisting of the following columns:
   * ```COMID```: Unique identifier 
    * ```fun_Btotal```: Change in functional LCBD of total beta diversity
    * ```fun_Brepl```: Change in functional LCBD of the replacement component
    * ```fun_Brich```: Change in functional LCBD of the richness difference component
    * ```phy_Btotal```: Change in phylogenetic LCBD of total beta diversity
    * ```phy_Brepl```: Change in phylogenetic LCBD of the replacement component
    * ```phy_Brich```: Change in phylogenetic LCBD of the richness difference component
    * ```tax_Btotal```: Change in taxonomic LCBD of total beta diversity
    * ```tax_Brepl```: Change in taxonomic LCBD of the replacement component
    * ```tax_Brich```: Change in taxonomic LCBD of the richness difference component
    * ```fun_Btotal_es```: Empirical effect size of change in ```fun_Btotal```
    * ```fun_Brepl_es```: Empirical effect size of change in ```fun_Brepl```
    * ```fun_Brich_es```: Empirical effect size of change in ```fun_Btotal```
    * ```phy_Btotal_es```: Empirical effect size of change in ```phy_Btotal```
    * ```phy_Brepl_es```: Empirical effect size of change in ```phy_Brepl```
    * ```phy_Brich_es```: Empirical effect size of change in ```phy_Brich```
* ```native_lcbd.rds```: Observed and ES values of taxonomic (no ES), functional, 
and phylogenetic local contributions to beta diversity (LCBD) for the native 
species pool. These data are used as explanatory variables in main analyses. 
Data frame consisting of the same structure as ```delta_lcbd.rds```, but LCBD 
values are for the native species, not a change over time.
* ```native_alpha.rds```: Observed and ES values of functional and phylogenetic 
richness for the native species pool. These data are used as explanatory 
variables in main analyses. Data frame consisting of the following columns:
    * ```COMID```: Unique identifier
    * ```fun_his_alpha```: Functional richness measured as the volume of kernel density hypervolumes.
    * ```phy_his_alpha```: Phylogenetic richness measured as the number of phylogenetic tree branches.
    * ```fun_his_alpha_es```: Empirical effect size of ```fun_his_alpha```
    * ```phy_his_alpha_es```: Empirical effect size of ```phy_his_alpha```

