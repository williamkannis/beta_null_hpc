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
library(purrr)

# Directories
dir <- commandArgs(trailingOnly = TRUE)
obs_dir <- file.path(dir,"obs_out")
null_dir <- file.path(dir,"null_out")
out_dir <- file.path(dir,"ses_outputs")
dir.create(out_dir,recursive = T)


# Load in obs and null data ----------------------------------------------------
obs_file <- list.files(obs_dir)
null_files <- list.files(null_dir)


obs_input <- readRDS(obs_file)
null_list <- lapply(
  null_files, 
  function(x) readRDS(file.path(null_dir,x))
)

## ADD IN CODE TO LEAVE MESSAGE IF DIRECTORIES HAVE FILES OTHER THAN THOSE NEEDED
# OR MAKE LIST FILE MORE SPECIFIC

# Combine null processing chunks into one list, with metrics outside, iterations
#inside
null_input <- purrr::transpose(null_list)


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
out <- purrr::transpose(out_list)


# Export -----------------------------------------------------------------------

out_name <-sub("_obs","",obs_file)
saveRDS(out,file.path(out_dir,out_name))

