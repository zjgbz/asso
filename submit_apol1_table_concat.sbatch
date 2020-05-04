#!/bin/bash

#SBATCH --job-name=apol1_table_concat
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00
#SBATCH --mem=5GB
#SBATCH --output=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.out
#SBATCH --error=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
#SBATCH --array=0-11

module purge
module load r/local_3.6.0

cohort_list=('GeneSTAR' 'GENOA' 'HyperGEN' 'WHI' 'ARIC' 'CARDIA' 'JHS' 'CHS' 'HCHS_SOL' 'FHS' 'MESA' 'DHS')
cohort=${cohort_list[${SLURM_ARRAY_TASK_ID}]}

status_name="apol1"
exp_int="expand_int"

Rscript apol1_table_concat.R ${status_name} ${cohort} ${exp_int} > ../${status_name}/${cohort}/result_table_apol1_${cohort}_${exp_int}.out 2>&1

exit
