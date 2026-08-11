#!/bin/bash

#SBATCH --job-name=phy_null_beta
#SBATCH --array=1-6
#SBATCH --cpus-per-task=37
#SBATCH --mem=45gb
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL

# Load in modules
module load anaconda3/2023.09-0
module load r/4.4.0
module load gdal/3.8.3
module load geos/3.12.1
module load proj/9.2.1
module load sqlite/3.43.2

# Run the task for each index in the job array
Rscript HPC/scripts_do_not_run/08_null_phy_beta.R HPC/beta_null_input_data/phy_null_input_list_${SLURM_ARRAY_TASK_ID}.rds
