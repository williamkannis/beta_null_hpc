#-------------------------------------------------------------------------------
#
#  Null model standardized effect size batch script
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 03/24/2026

# Description: Calculated standardized effects sizes, empirical p values, and
# diagnostic stats for a diversity metric defined by a shell script argument

## SET SHELL SCRIPT TO HAVE NUMBER OF COMPONENTD PLUS 1 IF LCBD = T

# House keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
# LOAD IN CUSTOM PACKAGE
library(dplyr)
library(purr)

# Directories
obs_dir <- file.path(dir,"obs_out")
null_dir <- file.path(dir,"null_out")
out_dir <- file.path(dir,"ses_outputs")
dir.create(out_dir,recursive = T)


# Load in obs and null data ----------------------------------------------------

obs <- readRDS(list.files(obs_dir))

# Combine null processing chunks into one list
null_files <- list.files(null_dir)
null_list <- lapply(
  null_files, 
  function(x) readRDS(file.path(null_dir,x))
)
null <- unlist(null_list,recursive = FALSE)


obs_list <- list(Beta = list(Btotal = 1,Brepl = 1),lcbd = 1)
null_list <- list(obs_list,obs_list,obs_list)

obs_input <- unlist(obs_list,recursive = FALSE)

null_list_t <- purrr::transpose(null_list)
null_input <- purrr::transpose(null_list_t[[1]])
## FIX DIV FUNCTION TO MAKE A CLEAN LIST ITH NO NESRING
if("lcbd" %in% names(null_list_t)) {
  null_input$lcbd <- null_list_t[[2]]
}


# Run SES function and export  -------------------------------------------------

out_list <- parallel::mclapply(1:length(obs_input), function(i){
  
  obs <- obs_input[[i]]
  null_list <- null_input[[i]]
  
  ses_fun(
    obs = obs,
    null_list = null_list,
    diag = T, 
    emp_padding =1
  )
}
)
names(out_list) <- names(obs_input)


# Export -----------------------------------------------------------------------
out_name <-sub("ses_input","ses_out",input_file)
saveRDS(out,file.path(out_dir,out_name))

