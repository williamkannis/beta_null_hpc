#-------------------------------------------------------------------------------
#
#  Taxonomic beta diversity observed cluster code
#
#-------------------------------------------------------------------------------

# Author 

# Created: 03/16/2026

# Description: Estimates taxonomic beta diversity and LCBD for for the observed
# values from the contemporary and native only species pools. For use on the 
# HPC cluster.

# House keeping ----------------------------------------------------------------
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
input_dir <- "HPC/beta_obs_input_data"
out_dir <- "HPC/obs_out"

# Load costume functions
source("diversity_batch_functions.R")

# Input data
time <- commandArgs(trailingOnly = TRUE)
input_file <- paste0(time,"_phy_obs_input_list.rds")
input_list <- readRDS(file.path(input_dir,input_file))
com <- input_list$com


# Estimate obs values ----------------------------------------------------------

# Run beta function
out <- tax_beta_batch(com = com)

# Export outputs
out_name <- paste0(time,"_tax_beta_obs",".rds")
saveRDS(out,file.path(out_dir,out_name))


# End of script

