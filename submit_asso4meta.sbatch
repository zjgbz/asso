#!/bin/bash

#SBATCH --job-name=asso4meta
##SBATCH -p general
#SBATCH -p bigmem
#SBATCH --qos bigmem_access
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=10GB
#SBATCH --output=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.out
#SBATCH --error=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.err
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=minzhi.hpc.status@gmail.com
#SBATCH --array=0-13

module purge
module load r/local_3.6.0

cohort_list=('HyperGEN' 'GeneSTAR' 'WHI' 'CHS' 'SAFS' 'DHS' 'GENOA' 'CARDIA' 'BioMe' 'COPDGene' 'ARIC' 'MESA' 'HCHS_SOL' 'JHS')
cohort=${cohort_list[${SLURM_ARRAY_TASK_ID}]}
# table_i=1

status="cn"

Rscript asso4meta.R ${status} ${cohort} table1 > ../${status}/${cohort}/asso4meta_table1_${cohort}_${SLURM_ARRAY_JOB_ID}.out 2>&1
Rscript asso4meta.R ${status} ${cohort} table2 > ../${status}/${cohort}/asso4meta_table2_${cohort}_${SLURM_ARRAY_JOB_ID}.out 2>&1
Rscript asso4meta.R ${status} ${cohort} table3 > ../${status}/${cohort}/asso4meta_table3_${cohort}_${SLURM_ARRAY_JOB_ID}.out 2>&1

exit
