#-------------------------------------------------------------------------------
#
#  Phylogenetic alpha diversity observed cluster code
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 03/16/2026

# Description: Estimates phylogenetic alpha diversity for the observed
# values from the contemporary and native only species pools. For use on  
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
out_dir <- "HPC/obs_out"

# Load costume functions
source("HPC/scripts_do_not_run/diversity_batch_functions.R")

# Input data
args <- commandArgs(trailingOnly = TRUE)
input_list <- readRDS(args)
com <- input_list$com
tree <- input_list$tree


# Estimate obs values ----------------------------------------------------------

# Run beta function
out <- BAT::alpha(comm = com,tree = tree)

# Export outputs
saveRDS(out,file.path(out_dir,"his_phy_alpha_obs.rds"))

