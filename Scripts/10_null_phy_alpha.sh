#!/bin/bash

#SBATCH --job-name=null_phy_alpha
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --time=30:00
#SBATCH --mail-type=ALL

# Load in modules
module load anaconda3/2023.09-0
module load r/4.4.0
module load gdal/3.8.3
module load geos/3.12.1
module load proj/9.2.1
module load sqlite/3.43.2

# Run the task for each index in the job array
Rscript HPC/scripts_hpc/10_null_phy_alpha.R HPC/alpha_null_input/his_phy_alpha_null_input_list.rds
