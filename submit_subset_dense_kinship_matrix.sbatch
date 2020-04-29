#!/bin/bash

#SBATCH --job-name=kinship_dense2sym
#SBATCH -p bigmem
#SBATCH --qos bigmem_access
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --mem=1200GB
#SBATCH --output=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A.out
#SBATCH --error=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com

module purge
module load r/local_3.5.0

Rscript kinship_dense2sym.R > result_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out 2>&1

exit
