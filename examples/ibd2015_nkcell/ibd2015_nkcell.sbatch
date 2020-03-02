#!/bin/bash

#SBATCH --job-name=ibd2015_nkcell
#SBATCH --output=/scratch/users/xiangzhu/dump/ibd2015_nkcell_%A_%a.out
#SBATCH --error=/scratch/users/xiangzhu/dump/ibd2015_nkcell_%A_%a.err
#SBATCH --array=1-125
#SBATCH --time=12:30:00
#SBATCH --nodes=1
#SBATCH --partition=hns,owners,normal,whwong
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xiangzhu.nku@gmail.com

module unload matlab
module load matlab/R2017b

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Copy the workhorse script to this sub-job
cp ibd2015_nkcell.m case_"$SLURM_ARRAY_TASK_ID"_run_analysis.m

# Add two lines in the sub-job script
echo "case_id = $SLURM_ARRAY_TASK_ID;" | cat - case_"$SLURM_ARRAY_TASK_ID"_run_analysis.m > temp_"$SLURM_ARRAY_TASK_ID"_run_analysis.m && mv temp_"$SLURM_ARRAY_TASK_ID"_run_analysis.m case_"$SLURM_ARRAY_TASK_ID"_run_analysis.m
echo "clear;" | cat - case_"$SLURM_ARRAY_TASK_ID"_run_analysis.m > temp_"$SLURM_ARRAY_TASK_ID"_run_analysis.m && mv temp_"$SLURM_ARRAY_TASK_ID"_run_analysis.m case_"$SLURM_ARRAY_TASK_ID"_run_analysis.m

# Run the sub-job script
matlab -nodisplay < case_"$SLURM_ARRAY_TASK_ID"_run_analysis.m

# Delete the sub-job script
rm case_"$SLURM_ARRAY_TASK_ID"_run_analysis.m
