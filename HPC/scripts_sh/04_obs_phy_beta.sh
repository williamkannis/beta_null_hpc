#!/bin/bash

#SBATCH --job-name=obs_phy_beta
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
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
Rscript HPC/scripts_hpc/04_obs_phy_beta.R HPC/beta_obs_input_data/${type}_phy_obs_input_list.rds
