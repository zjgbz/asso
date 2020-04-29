#!/bin/bash

#SBATCH --job-name=apol1_table
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=70GB
#SBATCH --output=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.out
#SBATCH --error=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
#SBATCH --array=1-4

module purge
module load r/local_3.6.0

cohort="APOL1"
pheno_i=2
exp_int="expand_int"

#Rscript formal_table${SLURM_ARRAY_TASK_ID}.R ${cohort} > ../cohort/${cohort}/asso_table${SLURM_ARRAY_TASK_ID}_${cohort}.out 2>&1
Rscript apol1_table.R ${cohort} ${SLURM_ARRAY_TASK_ID} ${pheno_i} ${exp_int} > ../cohort/${cohort}/asso_table3_${cohort}_${SLURM_ARRAY_TASK_ID}_${pheno_i}_${exp_int}.out 2>&1

exit
