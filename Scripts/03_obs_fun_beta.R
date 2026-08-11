#-------------------------------------------------------------------------------
#
#  Functional beta diversity observed cluster code
#
#-------------------------------------------------------------------------------

# Author: 

# Created: 03/16/2026

# Description: Estimates functional beta diversity and LCBD for for the observed
# values from the contemporary and native only speceis pools. For use on  
# HPC cluster.

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
out_dir <- "obs_out"

# Load costume functions
source("diversity_batch_functions.R")

# Input data
args <- commandArgs(trailingOnly = TRUE)
input_list <- readRDS(args)
com <- input_list$com
trait <- input_list$trait

# Output naming
time <- input_list$time


# Estimate obs values ----------------------------------------------------------

# Run beta function
out <- kernel_beta_batch(com = com,trait_hyp = trait)

# Export outputs
out_name <- paste0(time,"_fun_beta_obs",".rds")
saveRDS(out,file.path(out_dir,out_name))


# End of script

