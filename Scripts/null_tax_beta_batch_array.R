################################################################################
#
#  Tax beta diversity null models cluster code
#
################################################################################

# AUTHOR: William Annis

# CREATED: 02/14/2026

# DESCRIPTION: Batch estimates tax beta diversity and LCBD for a list of
# null model assemblages and exports files to directory. For use on the Palmetto
# cluster.

### House keeping  #############################################################
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
input_dir <- "beta_null_input_data"
out_dir <- "null_out"

# Load costume functions
source("diversity_batch_functions.R")

# Input data
cores <- 20
args <- commandArgs(trailingOnly = TRUE)
input_file <- paste0(args,"_tax_null_input_list.rds")
com_list <- readRDS(file.path(input_dir,input_file))


### Null model batch processing  ###############################################

# Run beta function
out_list <- mclapply(com_list,tax_beta_batch, mc.cores = cores)

# Check for scheduled core issue, do not export results if it exists
export_results <-all(sapply(out_list, function (x) length(x) == 2))

# Export outputs
if(export_results) {
  out_name <- paste0("mod_tax_beta_null.rds")
  saveRDS(out_list,file.path(out_dir,out_name))
} else {
  print("some scheduled cores did not return values. all jobs affected")
}


# End of script

