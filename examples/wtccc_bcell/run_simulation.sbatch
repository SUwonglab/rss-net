#!/bin/bash

#SBATCH --job-name=run_simulation
#SBATCH --output=/scratch/midway2/xiangzhu/tmp/dump/run_simulation_%A_%a.out
#SBATCH --error=/scratch/midway2/xiangzhu/tmp/dump/run_simulation_%A_%a.err
#SBATCH --array=0-1
#SBATCH --partition=broadwl
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32000 
#SBATCH --mail-type=ALL,ARRAY_TASKS 
#SBATCH --mail-user=xiangzhu.nku@gmail.com

# Load modules for RSS-NET
module unload matlab
module load matlab/2017b

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Copy the workhorse script to this sub-job
cp run_simulation.m case_"$SLURM_ARRAY_TASK_ID"_run_simulation.m

# Add two lines in the sub-job script
echo "trial_name = 'm$SLURM_ARRAY_TASK_ID';" | cat - case_"$SLURM_ARRAY_TASK_ID"_run_simulation.m > temp_"$SLURM_ARRAY_TASK_ID"_run_simulation.m && mv temp_"$SLURM_ARRAY_TASK_ID"_run_simulation.m case_"$SLURM_ARRAY_TASK_ID"_run_simulation.m
echo "clear;" | cat - case_"$SLURM_ARRAY_TASK_ID"_run_simulation.m > temp_"$SLURM_ARRAY_TASK_ID"_run_simulation.m && mv temp_"$SLURM_ARRAY_TASK_ID"_run_simulation.m case_"$SLURM_ARRAY_TASK_ID"_run_simulation.m

# Run the sub-job script
matlab -nodisplay < case_"$SLURM_ARRAY_TASK_ID"_run_simulation.m

# Delete the sub-job script
rm case_"$SLURM_ARRAY_TASK_ID"_run_simulation.m
