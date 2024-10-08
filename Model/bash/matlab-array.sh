#!/bin/bash

# Job Name
#SBATCH -J lca 

# Walltime requested
#SBATCH -t 144:00:00

# Provide index values (TASK IDs)
#SBATCH --array=0-120

#SBATCH --account=carney-ashenhav-condo

# Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID
##SBATCH -e job_files/arrayjob-%J-%a.err
#SBATCH -o ../jobfiles/arrayjob-%A-%a.out

# Controls the minimum/maximum number of nodes allocated to the job
#SBATCH -N 1

# single core
#SBATCH -c 1

#SBATCH --mail-type=ALL

#SBATCH --mail-user=xiamin_leng@brown.edu

# Use the $SLURM_ARRAY_TASK_ID variable to provide different inputs for each job
 

module load matlab/R2021a

echo "Running job array number 1: "$SLURM_ARRAY_TASK_ID

cd ../gridsearch/code

matlab-threaded -nodisplay -nojvm -r "LCA_effect_inhib_parallel($SLURM_ARRAY_JOB_ID,$SLURM_ARRAY_TASK_ID), exit"
