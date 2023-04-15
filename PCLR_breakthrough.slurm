#!/bin/bash

#SBATCH --job-name=PCLR
#SBATCH --mincpus=3
#SBATCH --mem-per-cpu=6GB
#SBATCH --time=2:30:00
#SBATCH --partition=standard
#SBATCH --account=guangyao
#SBATCH --array=1-1000

cd ---

module load R/4.2.2

Rscript PCLR_breakthrough.R ${SLURM_ARRAY_TASK_ID}
