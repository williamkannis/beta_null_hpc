#-------------------------------------------------------------------------------
#
#  Batch diversity HPC script
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 02/14/2026

# Description: Batch estimates observed or null iterations of user specified
# diversity metrics. Observed or null outputs are chosen based on shell script
# argument.

# House keeping  ---------------------------------------------------------------
rm(list = ls())


#Packages
# CAN REMOVE MOST OF THESE, AND JUST CALL IN THIS PACKAGE.
library(dplyr)
library(purrr)
library(tibble)
library(BAT)
library(adespatial)
library(ade4)
library(parallel)

# data
input <- readRDS()
type <- commandArgs(trailingOnly = TRUE)


# Estimate diversity  ----------------------------------------------------------

# Run beta function
input$type <- type
out_list <- do.call(null_iterations,input)


# Export outputs ---------------------------------------------------------------

# Export naming
dir <- input$dir  ## MAYBE CHANGW TO SHELL ARGUMNET
if (type == "obs") {
  out_dir <- file.path(dir,"obs_out")
  out_name <- paste0(facet,"_",metric,"_obs.rds")
}
if (type == "null") {
  
  # files are named by on consecutive iterations across all other nodes
  node_num <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  min_iter <- (node_num-1)*input$null_iter+1
  max_iter <- min_iter + input$null_iter -1
  out_dir <- file.path(dir,"null_out")
  out_name <- paste0(facet,"_",metric,"_null_",min_iter,"-",max_iter,".rds")
}
dir.create(out_dir,recursive = T)
saveRDS(out_list,file.path(out_dir,out_name))

# # WILL NEED TO FIGURE THE CHANGE VERSION
# time <- input_list$time
