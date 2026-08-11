#!/bin/bash

#SBATCH --job-name=tax_obs
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=00:05:00
#SBATCH --mail-type=ALL

# Load in modules
module load anaconda3/2023.09-0
module load r/4.4.0
module load gdal/3.8.3
module load geos/3.12.1
module load proj/9.2.1
module load sqlite/3.43.2

# Specify what dataset
types=("mod" "his")
type=${types[$SLURM_ARRAY_TASK_ID]}

# Run the task for each index in the job array
Rscript 02_obs_tax_beta.R HPC_data/beta_obs_input_data/${type}_fun_obs_input_list.rds
