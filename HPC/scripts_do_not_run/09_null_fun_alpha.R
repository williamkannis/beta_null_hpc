#-------------------------------------------------------------------------------
#
#  Functional alpha diversity null models cluster code
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 02/14/2026

# Description: Batch estimates functional alpha diversity for a list of
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
source("diversity_batch_functions.R")

# Input data
cores <- 2
args <- commandArgs(trailingOnly = TRUE)
input_list <- readRDS(args)
com_list <- input_list$com_list
trait <- input_list$trait

# Output naming
time <- input_list$time
iter <- names(trait_list)
min_iter <- iter[as.numeric(iter) == min(as.numeric(iter))]
max_iter <- iter[as.numeric(iter) == max(as.numeric(iter))]


# Null model batch processing --------------------------------------------------

# Run beta function
out_list <- mclapply(com_list,kernel_alpha_batch,trait = trait, mc.cores = cores)

out_name <- paste0("his_fun_alpha_null_",min_iter,"-",max_iter,".rds")
saveRDS(out_list,file.path(out_dir,out_name))

