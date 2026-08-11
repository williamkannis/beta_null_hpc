#-------------------------------------------------------------------------------
#
#  Functional beta diversity null models cluster code
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 02/14/2026

# Description: Batch estimates functional beta diversity and LCBD for a list of
# null model assemblages and exports files to directory. For use on HPC cluster


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

# Directories
out_dir <- "HPC/null_out"

# Load costume functions
source("HPC/scripts_do_not_run/diversity_batch_functions.R")

# Input data
cores <- 2
args <- commandArgs(trailingOnly = TRUE)
input_list <- readRDS(args)
com <- input_list$com
trait_list <- input_list$trait_list

# Output naming
time <- input_list$time
iter <- names(trait_list)
min_iter <- iter[as.numeric(iter) == min(as.numeric(iter))]
max_iter <- iter[as.numeric(iter) == max(as.numeric(iter))]


# Null model batch processing --------------------------------------------------

# Run beta function
out_list <- mclapply(trait_list,kernel_beta_batch,com = com, mc.cores = cores)

# Check for scheduled core issue, do not export results if it exists
export_results <-all(sapply(out_list, function (x) length(x) == 2))

# Export outputs
if(export_results) {
  out_name <- paste0(time,"_fun_beta_null_",min_iter,"-",max_iter,".rds")
  saveRDS(out_list,file.path(out_dir,out_name))
} else {
  print("some scheduled cores did return values. all jobs affected")
}

