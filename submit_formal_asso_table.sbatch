#!/bin/bash

#SBATCH --job-name=formal_asso_table
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=15GB
#SBATCH --output=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.out
#SBATCH --error=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
#SBATCH --array=0-14

module purge
module load r/local_3.6.0

cohort_list=('FHS' 'COPDGene' 'BioMe' 'ARIC' 'GeneSTAR' 'GenSalt' 'CARDIA' 'HyperGEN' 'DHS' 'HCHS_SOL' 'CHS' 'WHI' 'MESA' 'GENOA' 'JHS')
cohort=${cohort_list[${SLURM_ARRAY_TASK_ID}]}
# table_i=1

Rscript formal_asso_table.R ${cohort} table1 > ../cohort/${cohort}/asso_table1_${cohort}_${SLURM_ARRAY_JOB_ID}.out 2>&1
Rscript formal_asso_table.R ${cohort} table2 > ../cohort/${cohort}/asso_table2_${cohort}_${SLURM_ARRAY_JOB_ID}.out 2>&1
Rscript formal_asso_table.R ${cohort} table3 > ../cohort/${cohort}/asso_table3_${cohort}_${SLURM_ARRAY_JOB_ID}.out 2>&1

exit
