#!/bin/bash

#SBATCH --job-name=null_tax_beta
#SBATCH --cpus-per-task=20
#SBATCH --mem=100gb
#SBATCH --time=3:00:00
#SBATCH --mail-type=ALL

# Load in modules
module load anaconda3/2023.09-0
module load r/4.4.0
module load gdal/3.8.3
module load geos/3.12.1
module load proj/9.2.1
module load sqlite/3.43.2

# Run the task for each index in the job array
Rscript HPC/scripts_do_not_run/06_null_tax_beta.R
