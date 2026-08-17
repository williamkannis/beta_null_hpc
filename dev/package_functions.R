comm <- readRDS("dev/mod_com_diversity_input.rds")[1:10,-1]
trait <- readRDS("dev/trait_diversity_input.rds")
trait_gower <- BAT::gower(trait)
trait_hyp <- BAT::hyper.build(trait_gower,axes = 4)  # reduce using pcoa 
tree <- readRDS("dev/phylo_tree.rds")



# Input lit example  -----------------------------------------------------------
input <- list(
  null_iter = 10, # number of iterations per node
  null_cores = 10, # number of core per node
  label = "phy",  # tax, fun, phy, for file naming
  fun = ".dendrogram_beta",  # name of function to call in,
  fun_args = list(  # arugmenters for diversity function
    comm = as.matrix(comm),
    # trait = NULL,
    tree = tree,
    func = "sorenson",
    abund = F,
    comp = F
  ),
  lcbd = T,
  beta_comps = c("Btotal","Brepl","Brich"),
  algorithm = "taxa.labels"
)






#function testins

.kernel_alpha(comm = comm,trait = trait_hyp)
.kernel_beta(comm = comm,trait = trait_hyp,func = "jaccard",abund = F)
.hull_alpha(comm = comm,trait = trait_hyp)
.hull_beta(comm = comm,trait = trait_hyp,func = "jaccard")
.dendrogram_beta(comm,tree,func = "sorenson",abund = F, comp = F)
.dendrogram_alpha(comm,tree,func = "sorenson",abund = F, comp = F)

fun_args <- list(comm = comm,tree=tree,func = "sorenson",abund = F, comp = F)
do.call(.dendrogram_beta,fun_args)

fun <- ".dendrogram_beta"
fun <- ".dendrogram_alpha"
fun_args <- list(comm = comm,tree=tree,func = "sorenson",abund = F, comp = F)
fun_args <- list(comm = comm,tree=tree)

a <- null_iterations(type = "null",label = "phy",lcbd =T,
                beta_comps = c("Btotal","Brepl"),fun = fun,fun_args = fun_args,
                null_iter = 10,null_cores = 10)
purrr::transpose(a)

.div_fun(NULL,T,c("Btotal","Brepl","Brich"),fun,fun_args)

.div_fun("phy",T,c("Btotal","Brepl","Brich"),fun,fun_args)


# HPC script building  ---------------------------------------------------------


## WILL NEED ASCRIPT TO LOAD IN MODULES version numbers may cause problems
# SCRIPT TO INSTALL R PACKAGES ON CLUSTER

## WILL NEED INTIAL SCRIPT TO INSTALL ALL REQUIRED PACKAGE

#Function to estimate resources needed
estimate_hpc_resources <- function(){
  
  # Export
  make_shell_arg()
}

#function to create shell aurment tabl (optinal and used as helper)
make_shell_arg <- function(memory,cores,nodes,walltime) {
  data.frame(
    arg = c(
      "#SBATCH --array=1-",
      "#SBATCH --cpus-per-task=",
      "#SBATCH --mem=",  # WILL NEED TO ADD UNITS EVENTUALLY
      "#SBATCH --time=" # WILL NEED TO ADD UNITS EVENTUALLY
    ),
    value = c(
      memory,
      cores,
      nodes,
      walltime
    ),
    units = c(
      "",
      "",
      "gb",
      "")
  )
}
a <-make_shell_arg(18,2,999,40)
apply(a,1,paste,collapse = "")

# 1. Generate the individual lines
lines <- paste0(a$arg, a$value, a$units)
c <- c(
  "#!/bin/bash",
  "",
  lines
)
"Rscript HPC/scripts_do_not_run/diversity_calc.R obs"
"Rscript HPC/scripts_do_not_run/diversity_calc.R null"
cat(paste(c, collapse = "\n"))

# Collapse them into a single text block
slurm_block <- paste(c, collapse = "\n")

b <- readLines("HPC/scripts/01_obs_tax_beta.sh")

# View the raw string structure
cat(slurm_block)
# Function to make shell scripts
#  exports shells
#   obs
#     custom hpc args
#     arg = obs
#     custom dir
#   null
#   ses
#
#  exports data list
#   dir name
#   com
#   trait
#   phy
#   change?
#   alpha or beta
#   facet
#   function
#   function arguments
#   null model algorthm
args <-c(
  shell_dir,
  overwrite = F,
  estimate_hpc_args = T,
  hpc_args = NULL,
  change = F,   # T or F. estiamte change between two time steps. keep for later release
  method = "",  # taxonomic, dendrogram, hull, or kernel
  metric = "",  # alpha or beta
  label = "",  # could be tax, fun, phy, or anythin users specifes. used for file and column naming 
  algorithm = ""  # USE PICANTE NOTATION FOR THESE
  )
div_shell_builder <- function(
){
  
  # Input errors and warnings
  
  # create directory
  
  # check that all necassary function arguments are provided. MAYBE DONT KEEP 
  # THESE AS ..., MAYBE HAVE ALL VALABLE AND KEEP THEM AS NULL. WILL NEED TO EXPERIMENT
  
  # Create function name
  fun_name <- paste0(".",method,metric)
  
  # Estimate resource needs
  hpc_args <- estimate_hpc_resources()
  
  # create necassary shell scripts
  
  # add corresponding r scripts
  file.copy()
  
  # create input data
  ### Check that species are in correct order
  match <-trait_match(com,trait_hyp)
  com_match <- match[[1]]
  trait_match <- match[[2]]
  
  
  ###NEED TO HAVE ALL THE WARNINGS OF THE HPC SCRIPTS IN THIS FUNCTION AS WELL
  # THAT WAY USERS WILL NO THEIR INPUTS ARE INCORRECT
  
  
}

## SHell scripts

#1. observed

#2. null

#3. SES

### Local

## wrapper for beta deivrsity functions, allows for parrallel processing
parrallel_beta <- function(input,fun,...) {
  
  # named pair list
  pairlist
  
  # run function in parallel
  div_list <- parallel::mclapply(pairlist,function(x){
    
    # extract pair
    data <- input[x]
    
    # estimate diversity
    fun(data,...)
    
  })
  
  # Change into dist object
  
}

## local null functions?

## Effect size function

## algortihm functions

