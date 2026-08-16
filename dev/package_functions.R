comm <- readRDS("dev/mod_com_diversity_input.rds")[1:10,-1]
trait <- readRDS("dev/trait_diversity_input.rds")
trait_gower <- BAT::gower(trait)
trait_hyp <- BAT::hyper.build(trait_gower,axes = 4)  # reduce using pcoa 
tree <- readRDS("dev/phylo_tree.rds")

# Algo development  ------------------------------------------------------------
algorithm <- "frequency"
.comm_algorithm <- function(comm, algorithm) {
  picante_list <- eval(formals(picante::randomizeMatrix)$null.model)
  all_list <- c(picante_list,"taxa.labels")
  ## NEED A WARNING THOWN SOME WHERE FOR THIS TO STOP IF UNACCEPTABLE CHOICE IS
  # GIVEN
  if(algorithm %in% picante_list){
    null_comm <- picante::randomizeMatrix(comm,algorithm)
    return(null_comm)
  }
  
  if(algorithm == "taxa.labels") {
    null_comm <- comm
    colnames(null_comm) <- colnames(null_comm)[sample(ncol(null_comm))]
    return(null_comm)
  }
}

.trait_tree_swap(trait = NULL, tree = NULL) {
  if(!is.null(trait) & !is.null(tree)){
    stop("please provide either a trait matrix, or a tree. Not both.")
  }
  
  if(!is.null(trait)) {
    row.names(trait) <- row.names(trait)[sample(nrow(trait))]
    return(trait)
  }
  
  if(!is.null(tree)){
    picante::tipShuffle(tree)
  }
}

.comm_algorithm(comm,"taxa.labels")


# Input lit example  -----------------------------------------------------------
input <- list(
  null.iter = NULL, # number of iterations per node
  null.cores = NULL, # number of core per node
  label = "",  # tax, fun, phy, for file naming
  fun = "",  # name of function to call in,
  fun_agrs = list(  # arugmenters for diversity function
    comm = as.matrix(NULL),
    trait = as.matrix(NULL),
    tree = NULL,
    func = NULL,
    abund = NULL
  ),
  algorithm = ""
)

# Null model iteration function  -----------------------------------------------
null_iterations <- function(
    type,  # obs or null
    change = F,
    label = "",
    lcbd = NULL, # T or F
    beta_comps = NULL,  # all or vector of components matching BAT naming
    fun,  
    fun_args,
    algorithm = NULL,
    null.iter = NULL,
    null.cores = NULL
    ){
  
  if(!change){
    if(type == "obs"){
      obs <- .div_fun(
        label = label,
        lcbd = lcbd,
        beta_comps = beta_comps,
        fun = fun,
        fun_args = fun_args)
      return(obs)
    }
    
    if (type == "null"){
      
      # Set up null algorithm options
      args <- list(...)
      alg_input <- args[names(args) %in% c("comm","trait","tree")]
      alg_input <- c(alg_input,list(algorithm=algorithm))
      
      # estimate null iterations across processing cores
      null_iters <- parallel::mclapply(
        X = 1:null.iter,
        mc.cores = null.cores,
        FUN = function(i) {
          
          # Create null iteration
          null_input <- do.call(.null_algorithm,alg_input)
          null_args <- c(
            null_input,
            args[!names(args) %in% c("comm","trait","tree")]
          )
          
          # Estimate diversity
          .div_fun(
            label = label,
            lcbd = lcbd,
            beta_comps = beta_comps,
            fun = fun,
            fun_args = null_args
            )
        }
      )
    }
    return(null_iters)
  }
  # if(change){
  #   
  # }
}
# Helper functions  ------------------------------------------------------------

# estimate diversity using function of choice
.div_fun <- function(label = NULL,lcbd,beta_comps,fun,fun_args) {
  
  # Prepare custom labels if provided
  if(!is.null(label)) label <- paste0(label,"_")
  
  # Estimate diversity with specified method
  helper_fun <- get(fun)
  div <- do.call(helper_fun,fun_args)
  
  # Determine if alpha or beta
  metric = stringr::str_extract(fun, "(?<=_).*")
  
  if(metric == "alpha"){
    colnames(div) <- lapply(colnames(div),function(i) paste0(label,i))
    return(div)
  }
  
  # Select components and lcbd for beta diversity
  if (metric == "beta") {
    beta <- div[beta_comps]
    names(beta) <- lapply(names(beta),function(i) paste0(label,i))
    if (!lcbd) return(beta)
    
    # Create data frame with LCBD of all components
    lcbd <- purrr::map2(beta,names(beta),.lcbd_batch) %>% 
      purrr::reduce(left_join,by = join_by(COMID)) %>% 
      tibble::column_to_rownames("COMID")
    
    return(list(beta = beta,lcbd = lcbd))
  }
}


#Estimate LCBD and export as data.frame
.lcbd_batch <- function(D,name) {
  
  # check if euclidean
  sqrt <- ifelse(ade4::is.euclid(D) == T, F, T)
  
  # estimate LCBD
  lcbd <- adespatial::LCBD.comp(D,sqrt.D = sqrt)$LCBD
  
  # Create output dataframe
  out <- data.frame(COMID = labels(D))
  out[,name] <- lcbd
  row.names(out) <- out$COMID
  out
}

# Wrappers to easily call BAT functions
# FOR TESTING - TEST THAT VALUES EQUAL NON WRAPPER OUTPUTS
# AND THAT FUNCTIONS TAKE THE ARGUMENTS CORRECTLY
.kernel_alpha <- function(...){
  args <- list(...)
  kernel <- BAT::kernel.build(...)
  BAT::kernel.alpha(kernel)
}
.kernel_beta<- function(...){
  args <- list(...)
  kern_args <- args[names(args) != "func"]
  kernel <- do.call(BAT::kernel.build,kern_args)
  BAT::kernel.beta(kernel,func = args$func)
}

.hull_alpha <- function(...){
  hull <- BAT::hull.build(...)
  BAT::hull.alpha(hull)
}
.hull_beta <- function(...){
  args <- list(...)
  hull_args <- args[names(args) != "func"]
  hull <- do.call(BAT::hull.build,hull_args)
  BAT::hull.beta(hull,func = args$func)
}

.dendrogram_alpha <- function(...) BAT::alpha(...)
.dendrogram_beta <- function(...) BAT::beta(...)

.tax_alpha <-function(...) BAT::alpha(...)
.tax_beta <-function(...) BAT::beta(...)

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
                null.iter = 10,null.cores = 10)
.div_fun(NULL,T,c("Btotal","Brepl","Brich"),fun,fun_args)

.div_fun("phy",T,c("Btotal","Brepl","Brich"),fun,fun_args)

# HPC script building  ---------------------------------------------------------


## WILL NEED ASCRIPT TO LOAD IN MODULES version numbers may cause problems
# SCRIPT TO INSTALL R PACKAGES ON CLUSTER

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

