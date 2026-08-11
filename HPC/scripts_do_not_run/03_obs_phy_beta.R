#-------------------------------------------------------------------------------
#
#  Phylogenetic beta diversity observed cluster code
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 03/16/2026

# Description: Estimates phylogenetic beta diversity and LCBD for for the 
# observed values from the contemporary and native only species pools. For use 
# on the HPC cluster.


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
out_dir <- "HPC/obs_out"

# Load costume functions
source("diversity_batch_functions.R")

# Input data
args <- commandArgs(trailingOnly = TRUE)
input_list <- readRDS(args)
com <- input_list$com
tree <- input_list$tree

# Output naming
time <- input_list$time


# Estimate obs values ----------------------------------------------------------

# Run beta function
out <- phy_beta_batch(com = com,tree = tree)

# Export outputs
out_name <- paste0(time,"_phy_beta_obs",".rds")
saveRDS(out,file.path(out_dir,out_name))

