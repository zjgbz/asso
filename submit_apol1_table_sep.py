import time
import itertools
import subprocess
from subprocess import call
import os

def write_sbatch_script(Rscript_filename, pheno, job_name):
	with open('submit_%s.sbatch'%job_name, 'w') as out:
		out.write('#!/bin/bash\n')
		out.write('\n')
		out.write('#SBATCH --job-name=%s_%s\n'%(Rscript_filename, pheno))
		# out.write('#SBATCH -p general\n')
		out.write('#SBATCH -p bigmem\n')
		out.write('#SBATCH --qos bigmem_access')
		out.write('#SBATCH --nodes=1\n')
		out.write('#SBATCH --ntasks-per-node=1\n')
		out.write('#SBATCH --time=1-00:00:00\n')
		out.write('#SBATCH --mem=70GB\n')
		out.write('#SBATCH --output=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.out\n')
		out.write('#SBATCH --error=/proj/yunligrp/users/minzhi/output_error/asso/%x_%A_%a.err\n')
		out.write('#SBATCH --mail-type=END,ALL\n')
		out.write('#SBATCH --mail-user=minzhi.hpc.status@gmail.com\n')
		out.write('#SBATCH --array=1-10\n')
		out.write('\n')
		out.write('module purge\n')
		out.write('module load r/local_3.6.0\n')
		out.write('\n')
		out.write('status_name="apol1"\n')
		out.write('cohort="all_cohorts"\n')
		out.write('pheno_i=%s\n'%pheno)
		out.write('exp_int="expand_int"\n')
		out.write('rscript_filename="%s"\n'%Rscript_filename)
		out.write('\n')
		out.write('Rscript ${rscript_filename}.R ${status_name} ${cohort} ${SLURM_ARRAY_TASK_ID} ${pheno_i} ${exp_int} > ../${status_name}/${cohort}/result_${SLURM_ARRAY_TASK_ID}_${pheno_i}_${exp_int}.out 2>&1\n')
		out.write('\n')
		out.write('exit\n')
	return 1

def submit_job(Rscript_filename, pheno):
	job_name = "%s_%s"%(Rscript_filename, pheno)
	filename = "submit_%s.sbatch"%job_name
	write_sbatch_script(Rscript_filename, pheno, job_name)
	job_id = 1
	result = subprocess.run(['sbatch', '%s'%filename], stdout=subprocess.PIPE)
	job_id = result.stdout.decode('utf-8').rstrip().split(' ')[3]
	return job_id

def main():
	pheno_list = ["EGFRCKDEPI", "CKD"]
	Rscript_filename = "apol1_table_sep"

	job_id_list = []
	for pheno in pheno_list:
		tmp_job_id = submit_job(Rscript_filename, pheno)
		job_id_list.append(tmp_job_id)
	return job_id_list

if __name__ == '__main__':
	job_id_list = main()
	print(job_id_list)
