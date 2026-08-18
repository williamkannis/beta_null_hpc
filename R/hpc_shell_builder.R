args <-c(
  shell_dir,
  overwrite = F,
  estimate_hpc_args = T,
  hpc_args = NULL,
  label = NULL,
  algorithm = c("taxa.labels","frequency", "richness","independentswap", "trialswap"),
  iter = 999,
  change = F,   # T or F. estimate change between two time steps. keep for later release
  method = c("dendrogram","hull","kernel","taxonomic"),  
  metric = c("beta","alpha"),
  lcbd = NULL,
  beta_comps = NULL,
  ...
)
div_shell_builder <- function(
){
  
  # Input errors and warnings
  metric <- match.arg(metric)
  algorithm <- match.arg(algorithm)
  
  # create directory
  
  # check that all necassary function arguments are provided. MAYBE DONT KEEP 
  # THESE AS ..., MAYBE HAVE ALL VALABLE AND KEEP THEM AS NULL. WILL NEED TO EXPERIMENT
  
  # Create function name and check arguments
  fun <- .create_div_fun(method,metric,...)
  
  # Estimate resource needs
  hpc_args <- estimate_hpc_resources()
  
  # create necessary shell scripts
  .shell_construct()
  # add corresponding r scripts
  file.copy()
  
  # create input data
  ### Check that species are in correct order
  match <-trait_match(com,trait_hyp)
  com_match <- match[[1]]
  trait_match <- match[[2]]
  
  # Prepare hpc input file
  input <- list(
    null_iter = null_iter, # number of iterations per node
    null_cores = null_cores, # number of core per node
    label = label,  # tax, fun, phy, for file naming
    fun = fun,  # name of function to call in,
    fun_args = fun_args,  # arugmenters for diversity function
    lcbd = lcbd,
    beta_comps = beta_comps,
    algorithm = algorithm
  )
  saveRDS(input,file.path(dir,"hpc_input.rds"))
  
  
  ###NEED TO HAVE ALL THE WARNINGS OF THE HPC SCRIPTS IN THIS FUNCTION AS WELL
  # THAT WAY USERS WILL NO THEIR INPUTS ARE INCORRECT
  
  
}
fun_args
# NEED TO MAKE FUNCTION NOT ALLOW EXTRA ARGUMENTS, AND TO NOT ALLOW FUNCTIONS
# THAT DONT EXIST TO BE MADE
.create_div_fun(method="dendrogram",metric="alpha",comm = comm,tree=tree)
# Helper to create diveristy function name based on metric and method of choice.
# Also checks to ensure all proper arguents are supplied to the function
.create_div_fun <- function(method,metric,...){
  
  # Create function name
  fun_name <- paste0(".",method,"_",metric)
  
  # Prepare function arguments
  fun_args = list(...)  ## MAYBE MAKE COMM ITS OWN ARGUMENT AND HAVE FUNCTION TURN TO MATRIX
  
  # CHeck that all required arguments are given for specified function
  arg_names <- names(fun_args)
  fun_formals <- switch(
    fun_name,
    ".kernel_beta" = c(
      formals(BAT::kernel.beta),
      formals(BAT::hull.build),
    ),
    ".hull_beta" = c(
      formals(BAT::hull.beta),
      formals(BAT::hull.build)
    ),
    ".dendrogram_beta" = formals(BAT::beta),
    ".tax_beta" = formals(BAT::beta),
    ".kernel_alpha" = c(
      formals(BAT::kernel.alpha),
      formals(BAT::hull.build)
    ),
    ".hull_alpha" = c(
      formals(BAT::hull.alpha),
      formals(BAT::hull.build)
    ),
    ".dendrogram_alpha" = formals(BAT::alpha),
    ".tax_alpha" = formals(BAT::alpha)
  )
  
  required_args <- names(fun_formals[sapply(fun_formals, function(x) x == "")])
  if(!all(required_args %in% arg_names)) {
    miss_arg <- required_args[!required_args %in% arg_names]
    stop(paste0(
      "Please provide values for: ",
      do.call(paste, c(as.list(miss_arg), sep = ", "))
    )
    )
  }
  
  if(!all(names(fun_formals) %in% arg_names)) {
    formal_names <- names(fun_formals) 
    miss_arg <- fun_formals[!formal_names %in% arg_names]
    miss_names <- names(miss_arg)
    warning(paste0(
      "No values provided for: ",
      do.call(paste, c(as.list(miss_names), sep = ", ")),
      ". ",
      do.call(paste, c(as.list(miss_arg), sep = ", ")),
      " will be passed to function as defualt"
    )
    )
  }
  fun_name
}
