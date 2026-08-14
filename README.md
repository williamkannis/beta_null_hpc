# Beta diversity change null modelling workflow for high performance computer (HPC) clusters

# Anonymous peer review version

This repository contains a R and Slurm workflow for calculating observed and 
null model standardized alpha diversity, beta diversity, and local contribution 
to beta diversity (LCBD) changes. Workflow estimates values in the taxonomic, 
functional, and phylogenetic dimensions.

Estimating multiple null iterations of multidimensional diversity is 
computationally demanding and can have large wall times, this is especially the 
case with beta diversity. As such, This workflow is designed to run using 
R on local machines  with the more expensive calculations ran using high 
performance computing (HPC) clusters via the slurm interface and shell scripts.

The workflow was developed for the two following analyses:

[Null models reveal differing drivers of multidimensional beta diversity change in invaded metacommunities](BLANK)

[Invasion syndromes based on shared changes in multidimensional beta diversity](BLANK)

Refer to these studies for detailed methodology and justification.

# Citation and contact information
If you use or adapt this workflow, please cite:

BLANK. Null models reveal differing drivers of multidimensional beta diversity 
change in invaded metacommunities. in review

BLANK. Invasion syndromes based on shared changes in multidimensional beta diversity.
in review

[Workflow citation / DOI]


For questions about this analysis, please contact:

```bash
Name: William K. Annis

Email: wannis@fsu.edu, williamkannis@gmail.com

OrcID: 0009-0003-3541-8503
```

# Purpose
The goal of this repository is 1) provide code used in the analyses from the
associated manuscripts, and 2) serve as general workflow to guide future null
model analyses

We are not able to directly share the community data form the manuscripts required
to replicate the manuscript's beta diversity metrics. As such, this repository 
serves as a guideline for performing beta diversity change null model analyses 
on HPC clusters. To reproduce the results of the associated
manuscripts using provided intermediate data, see the following repositories:

[Null models reveal differing drivers of multidimensional beta diversity change in invaded metacommunities](BLANK)

[Invasion syndromes based on shared changes in multidimensional beta diversity](BLANK)


## Workflow overview

```text
                            Input data
                                │
                                ▼
                        Prepare input data
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
                    Summarize null model results
                                │
                                ▼
                        Final effect-size data
```
This workflow estimates observed and null model effect sizes of multidimensional
alpha, beta, and LCBD diversity. Users prepare 
[diversity input data](#input-data) on their local machine and run the 
observed and null iteration calculations on HPC clusters. Observed and null
iterations are used to calculated standardized effect sizes on the HPC cluster.
Users can then download the effect size data back to their local machine to 
summarize or conduct analyses on these data.

While this workflow is tailored towards the manuscripts' data structure and 
methodology, we offer recommendations to generalize the workflow 
[here](#generalizing-this-workflow).


## Workflow null model end products
<table>
  <thead>
    <tr>
      <th rowspan="2">Diversity Facet</th>
      <th colspan="3">Native</th>
      <th colspan="3">Contemporary</th>
      <th colspan="3">Contemporary - native</th>
    </tr>
    <tr>
      <th>Alpha</th>
      <th>Beta</th>
      <th>LCBD</th>
      <th>Alpha</th>
      <th>Beta</th>
      <th>LCBD</th>
      <th>Alpha</th>
      <th>Beta</th>
      <th>LCBD</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Taxonomic</td>
      <td>—</td>
      <td>—</td>
      <td>—</td>
      <td>—</td>
      <td>—</td>
      <td>—</td>
      <td>—</td>
      <td>✓</td>
      <td>✓</td>
    </tr>
    <tr>
      <td>Functional</td>
      <td>✓</td>
      <td>✓</td>
      <td>✓</td>
      <td>—</td>
      <td>✓</td>
      <td>✓</td>
      <td>—</td>
      <td>✓</td>
      <td>✓</td>
    </tr>
    <tr>
      <td>Phylogenetic</td>
      <td>✓</td>
      <td>✓</td>
      <td>✓</td>
      <td>—</td>
      <td>✓</td>
      <td>✓</td>
      <td>—</td>
      <td>✓</td>
      <td>✓</td>
    </tr>
  </tbody>
</table>

The two associated manuscripts used a pseudo-historical approach to examine 
patterns and drivers of multidimensional beta diversity change. As such it was
necessary to estimate functional and phylogenetic pairwise beta diversity and LCBD 
at two species pools representing current (contemporary) and historic 
(native-only) communities, and estimate the change in diversity between these 
species pools.

This workflow also produces an optional null model output for changes in taxonomic 
diversity and LCBD, using an [algorithm](#taxonomic-beta-diveristy-change) that 
shuffles that occurrences of nonnative species. This null model does not 
produce separate null iterations for native
and contemporary communities, and is only used for change values. The taxonomic
null model is only suited for pseudo-historical approaches, where as the functional
and phylogenetic aspects of the workflow can work with true historic data.

Additionally, native functional and phylogenetic alpha and beta diversity were 
estimated as potential drivers of beta diversity change. Observed diversity 
values are produced for every diversity metric of interest.


See [Output data structure](#output-data-strucutre) for more information on
data products.

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
│   ├── mod_com_diversity_input.rds*
│   ├── his_com_diversity_input.rds*
│   ├── trait_diversity_input.rds*
│   └── phylo_tree.rds*
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

(*) Users will need to add these files
```
Place community data for both species pools (contemporary - "mod" and native - 
"his"), trait data, and phylogenetic trees in ```Diversity Input Data``` with the 
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

Scripts needing to be ran by the user can be found in ```hpc/scripts``` 
directory and are number based on order of workflow. Scripts found
in ```hpc/scripts_do_not_run``` are either scripts ran on the HPC cluster using 
shell scripts, or contain functions used by other scripts. **DO NOT** directly run 
scripts in ```hpc/scripts_do_not_run```. See [here](#generalizing-this-workflow) 
for tips on generalizing workflow to other data.


### 1. Prepare input data
**Perform on local machine**

**Script:** ```00_null_input_creation.R```

<ins>Purpose:</ins> Prepares input data for the calculation of observed and null
diversity values on high performance computer clusters. By default, 999 iteration
of random community matrices (regionally constrained and taxonomic beta null 
models), random functional-spaces/phylogenetic-trees (taxa-swap null models).
See the [Null model methods](#null-model-methods) section for information on 
null model algorithms used as default. Input data are prepared in chunks of 
random iterations to run on multiple PC nodes simultaneously. Each chunk 
is ran on one node, with parallel processing occurring within node. Workflow
contains the number of chunks and computational resources used in the associated
manuscripts. For more information on HPC resources and recommended data chunking
for general use of the workflow, see [HPC resources](#1-hpc-resources) section.

<ins>Outputs:</ins> Lists of input data for each metric and processing chunk.
For all metrics but taxonomic, list contains community and 
functional/phylogenetic data. Depending on null model, either the community or
functional/phylogenetic data are a list of 999 random iterations.


### 2. Upload entire ```HPC``` directory to high performance cluster storage.
For all steps utilizing HPC clusters, users will run shell scripts that will
run their respective R script using specified HPC resources. Below we only
list the shell scripts, but R scripts can be found in ```scripts_do_not_run``` and
have the same name as their shell scripts. 


### 3 .Calculate observed diversity values
**Perform on HPC**

**Scripts:** ```01_obs_tax_beta.sh```, ```02_obs_fun_beta.sh```
, ```03_obs_phy_beta.sh```, ```04_obs_fun_alpha.sh``` , ```05_obs_phy_alpha.sh```

<ins>Purpose:</ins> Calculates the observed beta diversity for the contemporary 
and native only species pools, and the observed alpha diversity of the native 
only species pool. Uses one high performance computer nodes for each time step. 

<ins>Outputs:</ins>

* [Observed beta diversity data](#observed-beta-diversity-data)
* [Observed alpha diversity data](#observed-alpha-diversity-data)

### 4. Calculate null diversity values
**Perform on HPC**

**Scripts:** ```06_null_tax_beta.sh```, ```07_null_fun_beta.sh```
, ```08_null_phy_beta.sh```, ```09_null_fun_alpha.sh```, ```10_null_phy_alpha.sh```

<ins>Purpose:</ins> Calculates the null iterations for beta diversity in the 
contemporary and native only species pools, and null iterations for alpha 
diversity in the native only species pool. 

<ins>Outputs:</ins> [Null iterations](#null-iterations)


### 5. Prepare effect size input data
**Perform on HPC**

**Scripts:** ```11_tax_beta_null_model_prep.sh```, ```12_beta_null_model_prep.sh```
, ```13_alpha_null_model_prep.sh```

<ins>Purpose:</ins> Calculates difference in LCBD between contemporary and 
native pools (delta) for the observed values, and for each null iteration. 
Consolidates outputs into single files, with separate delta, native, and 
contemporary species pool values.

<ins>Outputs:</ins> [Intermediate files used for effect size calculations](#formatted-ses-inputs)


### 6. Calculate effect sizes
**Perform on HPC**

**Script:**```14_batch_ses.sh```

<ins>Purpose:</ins> Estimates null model effect sizes of each metric, as well
as null model diagnostics. See [here](#effect-size-calculations) 
for more information on effect size calculations. Runs on 1 node per diversity 
metrics of interest  (i.e., number of items in ```ses_inputs``` directory).

<ins>Outputs:</ins> [Summarized null model outputs](#summarized-null-model-outputs)

### 7. Download the entire ```HPC``` directory to local machine.


### 8. Summarize null model results
**Perform on local machine**

**Scripts:** ```15_ses_comp.R```

<ins>Purpose:</ins> Compiles and formats the resulting SES, ES, and diagnostic 
stats across files. Also visualizes normality diagnostics to allow users to 
decide between SES or ES values for further analyses.

<ins>Outputs:</ins> 

* [```delta_lcbd.rds```](#effect-size-dataframes), 
* [```native_lcbd.rds```](#effect-size-dataframes), 
* [```native_alpha.rds```](#effect-size-dataframes)


### Helper function scripts 
**Directory:**```scripts_do_not_run```

* ```null_model_algorithms.R```: Contains algorithms to randomize community, 
trait, or phylogenetic data for null model analysis
* ```diversity_batch_functions.R```: Contains functions that estimate multiple 
iterations of diversity metrics using parallel computation.
* ```null_model_effect_size_function.R```: Function to summarize null 
distributions and calculates standardize effect sizes in the traditional z-score 
method (SES), empirical p-values, and p-value based effect sizes (ES), and 
reports optional diagnostic metrics.


## Output data strucutre
See below for detailed descriptions of intermidate and final data products from
this workflow.

### Observed beta diversity data 
**Directory:**```HPC/obs_out```

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
**Directory:**```HPC/obs_out```

**Files:**```his_fun_alpha_obs.rds```, ```his_phy_alpha_obs```

Files contain named numeric objects with observed functional richness 
(volume of kernel density hypervolume) or phylogenetic richness (number of 
branches in phylogenetic tree). 

### Null iterations 
**Directory:**```HPC/null_out```

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
**Directory:**```HPC/ses_input```

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
* ```$null_list```: list of null iterations
   
### Summarized null model outputs 
**Directory:**```HPC/ses_output```

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
**Directory:**```Diversity Output Data```

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

## Null model methods
The functions in ```null_model_algorithms.R``` represent the null model 
algorithms used in the manuscripts associated with this workflow. Below we
briefly describe the algorithms used in the workflow

### Taxonomic beta diveristy change
The purpose of this null model was to determine if observed patterns of beta
diversity change reflected non-random patterns of nonnative species occurrences.
To do so we used a homogenization null model developed by Leprieur et al. (2007), 
for the use in studies that use contemporary and native species pools
as proxies for true historical data. This null model maintains an equiprobable 
total of columns by randomizing the spatial distribution of nonnative species, 
maintaining their occurrence frequency, but allowing communities to be 
unconstrained in the number of nonnative species they receive. As our the 
manuscripts associated with this workflow covers large spatial extents with 
biogeographical barriers, we constrained the model to only allow nonnative 
species to be assigned to sites within the regional species pools (hydrological 
regions, HUC2s) in which they have been introduced. 
See [here](https://doi.org/10.1111/j.1472-4642.2007.00409.x) 
for more information and arguments for the ecological validity of this null model 

### Phylogenetic and functional beta diversity
Functional and phylogenetic diversity metrics are often inherently linked to 
underlying taxonomic structure and require null model standardization to 
determine the effects of species traits/phylogenies in diversity change. 
When disentangling the effects of taxonomic structure on functional and 
phylogenetic beta diversity metrics, it is important to maintain taxonomic 
beta diversity. Therefore, we used a taxon-swap algorithm according to the 
entire study extent as a species pool to preserve the observed taxonomic beta 
diversity. The taxon-swap algorithm swaps the species labels on trait matrices 
and phylogenetic trees to randomize traits/phylogeny while maintaining species 
richness, occurrence frequencies, and co-occurrences 

### Phylogenetic and functional alpha diversity
For native alpha diversity, we used a regionally-constrained taxon-swap model to 
reflect historical dispersal limitations caused by biogeographic differences 
among watershed boundaries. Here, we randomized community matrices, shuffling 
species within regional species pools, which maintained species richness, 
occurrence frequencies, and co-occurrences within regional species pools. 



## Effect size calculations


Using ```null_model_effect_size_function.R```, this workflows estimates null 
model effect sizes in two different ways:


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

Where ```obs``` refers to observed diversity values, ```null``` refers to
null distribution, ```sd``` refers to standard devation, ```p``` refers to the
emprial p-value, ```r``` refers to the number of null iterations more extreme
than the observed values, and ```n``` refers the the number of null iterations.

This functional also summarizes the null distributions, and reports optional 
noramility diagnostics used to select between the two effect size methods. 
Asymmetrical null distributions should be assessed using empirical p-value based 
effect sizes rather than z-score based SES. See 
[Botta-Dukát (2018)](https://doi.org/10.1556/168.2018.19.1.8) for more 
information on selecting SES or p-value based ES.

## Generalizing this workflow
The workflow as provided is designed to produce the beta diversity metrics used 
in the two associated manuscripts, but offer an overall general framework. User
will understandably need to make some adjustments to the workflow to match their
data and goals. The diversity output data structure and most of the workflow 
will remain the same for most use cases. However, there are certain cases that 
will result in changes to the workflow:

1. Users using this workflow with their own data, regardless of structure, 
will likely require different HPC resources (e.g. nodes, cores, memory).
2. Users may wish to use different null model algorithms.
3. Users may wish to only a subset of the possible diversity metrics
4. Users may either want to use a single time step, not estimate change, or use
more than 2 time steps.

Below, we outline how users will need to adapt workflow to accommodate these 
changes

### 1. HPC resources
Estimating multiple null iterations of multidimensional diversity is 
computationally demanding and can have large wall times, this is especially the 
case with beta diversity. To reduce this wall time, our workflow runs multiple 
null iterations simultaneously within and among HPC nodes.

Input data are prepared in chunks of random iterations to run on multiple HPC nodes
simultaneously, with each chunk being ran on one node. Each node will be able
to run multiple null model iterations simultaneously across a user specified
amount of cores. By separating iterations across nodes, users can use more 
computer cores simultaneously, reduce memory pressure, reduce HPC queue times,
and more efficiency use HPC resources.

In the workflow, the default number of chunks and iterations are based on the 
memory requirement and needs of the associated manuscripts. Consequently, the number of nodes and
cores in the workflow are also based on these requirements. For example, the 1998 
iterations (999 iterations each pool) of functional beta diversity we divided 
into 999 chunks with 2 iteration per node due to the high memory requirements of 
kernel density based functional metrics. Conversely, the 1998 iteration of 
phylogenetic beta diversity were divided into 6 chunks with 333 iterations per 
chunk

Users should determine chunk sizes based their own computational needs, and
can choose the number of chunks per species pool, and number of iterations by 
altering the following code in ```00_null_input_creation.R```:

```bash
max_iter  <- 999
```

```bash
# How many cores per nodes to get 999 total
sapply(c(1, 3, 9, 27, 37, 111),function(x) 999/x)

# Create splits based on number of nodes needed
n_groups_t <- 1     # taxonomic
n_groups_f <- 500  # functions
n_groups_p <- 3    # phylogenetic
```

Users will then need to alter the shell scripts in the 
[observed](#3-calculate-observed-diversity-values) 
and [null iteration](#4-calculate-null-diversity-values) 
workflow steps to reflect the number of chunks, required memory, number cores, 
and wall time for each metric. See below for an example:

```bash
# Process 6 chunks using 6 HPC nodes 
#SBATCH --array=1-6  

# Use 37 cpu cores per node (i.e., 37 parralel operations per chunk)
#SBATCH --cpus-per-task=37

# Use 45gb of ram per node
#SBATCH --mem=45gb

# 6 hours of walltime
#SBATCH --time=6:00:00
```
<ins>TIP:</ins> We recommend that users experiment with memory and CPU 
requirements with smaller number of null iterations before running full job.

<ins>TIP:</ins> Do not set cpu cores to value higher than number of iterations in 
each chunk, as this will waste resources. For example, if 2 iterations per chunk,
the maximum allowable number of cores is ```SBATCH --cpus-per-task=2 ```.

### 2. Null model algorithms choices
The current workflow is limited in null model [algorithms](#null-model-methods) 
to those used in the associated manuscripts. While these algorithms fit many
broad cases such as functional and phylogenetic beta diversity, there are some
limitations with the alpha diversity and taxonomic beta diversity null models.
For example, users measuring taxonomic beta diversity change with true 
historical data could benefit from an independent swap algorithm. Users should 
research the best null  model algorithm for their usage and modify ```00_null_input_creation.R``` 
and ```null_model_algorithms.R``` accordingly. Future releases will contain more
null model algorithms.

### 2. Diversity output choices
The current workflow estimates taxonomic, functional, and phylogenetic alpha,
diversity, beta diversity, and LCBD. If users only want to measure certain
dimensions (e.g., phylogenetic only), or a certain metric (e.g., no alpha
diversity), then users should ignore the corresponding sections 
of ```00_null_input_creation.R``` that produce input data for those metrics.

Additionally, users should only run the shell scripts in ```HPC/scripts``` that
correspond to the dimensions and metrics of interest. See table below for script 
naming codes:

```bash
tax = taxonomic
fun = functional
phy = phylogenetic

obs = observed
null = null iterations
```
Ignoring certain metrics will not require additional changes to the workflow.
Users will just need to adjust ```#SBATCH --array=1-26``` 
in ```14_ses_batch.sh``` to reflect the number of diversity
metrics in ```HPC/ses_inputs```.

### 3. Species pools and change
This workflow can easily accommodate single time steps of diversity, or to 
not estimate beta diversity change. This will only require the modification
of ```00_null_input_creation.R```, and ```11_tax_beta_null_model_prep.R``` 
or ```12_beta_null_model_prep.R``` scripts

Users wanting only to estimate one time step should ignore all code in the above 
scripts that creates objects prefixed by "his_" (native pool). Additionally,
users should avoid code that creates objects prefixed by "d_" if they do not
want to create beta diversity change metrics.

Adding multiple time steps with no change will require users to create new objects
for each time step in the above scripts and ignoring the "d_" object codes. Users
should use the existing code as a guide.

Adding multiple time steps with change between each will require heavy modification
of the ```11_tax_beta_null_model_prep.R``` and ```12_beta_null_model_prep.R``` 
scripts to include functions to estimate multiple change comparisons.
