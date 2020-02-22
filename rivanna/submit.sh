#!/bin/bash
#SBATCH -A janeslab
#SBATCH -p standard
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -o ./logs/%A_%a.out

module load gcc R

Rscript ./leaveOneOut.R $SLURM_ARRAY_TASK_ID
