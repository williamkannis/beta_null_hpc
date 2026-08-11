#-------------------------------------------------------------------------------
#
#  Alpha Null Model Preparation
#
#-------------------------------------------------------------------------------

# Author: 

# Created: 03/24/2026

# Description: Prepares null and observed values for each diversity metric into
# format that can be easily summarized. Null model files are compiled into one,
# file. This script is for use on a HPC cluster and import functional or 
# taxonomic diversity values based on shell script argument.

# House keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(dplyr)
library(purrr)
library(parallel)

# Directories
obs_dir <- "obs_out"
null_dir <- "null_out"
out_dir <- "ses_inputs"

# Load in data  ----------------------------------------------------------------

# Select dimension to import
arg <- as.numeric(commandArgs(trailingOnly = TRUE))  # shell script argument
dim_list <- c("fun","phy")
dim <- dim_list[arg]
his_dim <- paste0("his_",dim,"_alpha")

# Load in observed data
his_obs <- readRDS(file.path(obs_dir,paste0(his_dim,"_obs.rds")))

# Extract null files
his_null_files <- list.files(null_dir,(his_dim))

# Load in null files
his_null_list <- lapply(his_null_files, function(x) readRDS(file.path(null_dir,x)))

# Combine null processing chunks into one list
his_null <- unlist(his_null_list,recursive = FALSE)

# Export -----------------------------------------------------------------------
out_list <- list(obs=his_obs,null_list=his_null)
out_name <- paste0(dim,"_his_alpha_ses_input.rds")
saveRDS(out_list,file.path(out_dir,out_name))

