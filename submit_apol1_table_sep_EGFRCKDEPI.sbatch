#!/bin/bash

#SBATCH --job-name=apol1_table_sep_EGFRCKDEPI
#SBATCH -p bigmem
#SBATCH --qos bigmem_access#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=70GB
#SBATCH --output=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.out
#SBATCH --error=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
#SBATCH --array=1-10

module purge
module load r/local_3.6.0

status_name="apol1"
cohort="all_cohorts"
pheno_i=EGFRCKDEPI
exp_int="expand_int"
rscript_filename="apol1_table_sep"

Rscript ${rscript_filename}.R ${status_name} ${cohort} ${SLURM_ARRAY_TASK_ID} ${pheno_i} ${exp_int} > ../${status_name}/${cohort}/result_${SLURM_ARRAY_TASK_ID}_${pheno_i}_${exp_int}.out 2>&1

exit
