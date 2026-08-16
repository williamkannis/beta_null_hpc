#-------------------------------------------------------------------------------
#
#  Batch diveristy HPC script
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 02/14/2026

# Description: Batch estimates observed or null iterations of user specified
# diversity metrics. Observed or null outputs are choosen based on shell script
# argument.

# House keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(dplyr)
library(purrr)
library(tibble)
library(BAT)
library(adespatial)
library(ade4)
library(parallel)


# Load custom functions
source("HPC/scripts_do_not_run/diversity_batch_functions.R")

# data
input <- readRDS()
type <- commandArgs(trailingOnly = TRUE)
input_list <- readRDS(args)
com <- input_list$com
trait_list <- input_list$trait_list

# Output naming
time <- input_list$time
iter <- names(trait_list)
min_iter <- iter[as.numeric(iter) == min(as.numeric(iter))]
max_iter <- iter[as.numeric(iter) == max(as.numeric(iter))]


# Estimate diversity  ----------------------------------------------------------

# Run beta function
out_list <- null_iterations()


# Export outputs ---------------------------------------------------------------

# Export naming
dir <- input$dir
if (type == "obs") {
  out_dir <- file.path(dir,"obs_out")
  out_name <- paste0(facet,"_",metric,"_obs.rds")
}
if (type == "null") {
  node_num
  min_iter <- (node_num-1)*input$null_iter+1
  max_iter <- min_iter + input$null_iter -1
  out_dir <- file.path(dir,"null_out")
  out_name <- paste0(facet,"_",metric,"_null_",min_iter,"-",max_iter,".rds")
}
saveRDS(out_list,file.path(out_dir,out_name))


# # NEED TO THINK ABOUT CHANGE VERSION
# time <- input_list$time
# out_name <- paste0(time,"_fun_beta_null_",min_iter,"-",max_iter,".rds")
# 
# # DO YOU STILL WANT THIS CHECK?
# # Check for scheduled core issue, do not export results if it exists
# export_results <-all(sapply(out_list, function (x) length(x) == 2))
# 
# # Export naming
# # Export outputs
# if(export_results) {
#   out_name <- paste0(time,"_fun_beta_null_",min_iter,"-",max_iter,".rds")
#   saveRDS(out_list,file.path(out_dir,out_name))
# } else {
#   print("some scheduled cores did return values. all jobs affected")
# }
# 
