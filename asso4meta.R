library(gaston)
library(feather)
library(dict)

invnorm <- function(expression){
  return(qnorm((rank(expression, na.last="keep") - 0.5)/sum(!is.na(expression))))
}

load_kinship <- function(load_dir, filename_prefix) {
	kinship_dir_filename_prefix = file.path(load_dir, filename_prefix)
	kinship_dir_filename_tsv = sprintf("%s.tsv", kinship_dir_filename_prefix)
	kinship_dir_filename_feather = sprintf("%s.feather", kinship_dir_filename_prefix)
	if (file.exists(kinship_dir_filename_tsv)) {
		kinship = read.csv(kinship_dir_filename_tsv, sep='\t', header=TRUE, row.names=1)
	} else if (file.exists(kinship_dir_filename_feather)) {
		kinship = read_feather(kinship_dir_filename_feather)
		kinship = as.matrix(apply(kinship, 2, as.numeric))
	}
	return (kinship)
}

args = commandArgs(trailingOnly=TRUE)
status = args[1]
cohort = args[2]
table_i = args[3]

biclass_list = list("CKD", "microcytosis", "anemia")

table_dict = dict()
table_dict[["table1"]] = list("cn-binary_cn0", "cn-binary_cn1", "cn-binary_cn3", "rs334", "rs33930165")
table_dict[["table2"]] = list("rs334_del0", "rs334_del1", "rs334_del2", "del-rs334")
table_dict[["table3"]] = list("rs11248850", "rs11248850_del0", "rs11248850_del12", "rs11248850_cn2", "rs11248850_cn34", "del-rs11248850")
var_list = table_dict[[table_i]]

pheno_list = list("hemoglobin_mcnc_bld_1", "hematocrit_vfr_bld_1", "rbc_ncnc_bld_1", "DDIMER",
                  "mcv_entvol_rbc_1", "mch_entmass_rbc_1", "mchc_mcnc_rbc_1", "rdw_ratio_rbc_1",
                  "neutrophil_ncnc_bld_1", "lymphocyte_ncnc_bld_1", "basophil_ncnc_bld_1",
                  "eosinophil_ncnc_bld_1", "monocyte_ncnc_bld_1", "wbc_ncnc_bld_1", "lnHBA1C", 
                  "pmv_entvol_bld_1", "platelet_ncnc_bld_1", "EGFRCKDEPI", "CKD", "microcytosis", "anemia")

header_list = list()
for (var_i in var_list) {
	if (grepl("-", var_i) && !grepl("cn-binary", var_i)) {
		header_list = append(header_list, sprintf("N_%s", var_i))
		header_list = append(header_list, sprintf("p-val_%s", var_i))
	} else {
		header_list = append(header_list, sprintf("N_%s", var_i))
		header_list = append(header_list, sprintf("beta_%s", var_i))
		header_list = append(header_list, sprintf("SE_%s", var_i))
		header_list = append(header_list, sprintf("p-val_%s", var_i))
	}
}

row_num = length(pheno_list)
col_num = length(header_list)
table_df = data.frame(matrix(ncol = col_num, nrow = row_num))
colnames(table_df) = header_list
row.names(table_df) = pheno_list

for (var_i in var_list) {
	for (pheno_i in pheno_list) {
		print(var_i)
		print(pheno_i)
		load_dir = file.path('..', status, cohort, table_i, var_i)
		y_dir_filename = sprintf('../%s/%s/%s/%s/common_pheno_%s_%s_%s.tsv', status, cohort, table_i, var_i, table_i, var_i, pheno_i)
		if (!(file.exists(y_dir_filename))) {
			table_df = table_df[!(row.names(table_df) %in% c(pheno_i)), ]
			next
		} else if (file.exists(y_dir_filename)) {
			y = read.csv(sprintf('../%s/%s/%s/%s/common_pheno_%s_%s_%s.tsv', status, cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
			var = read.csv(sprintf('../%s/%s/%s/%s/common_var_%s_%s_%s.tsv', status, cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
			age_sex = read.csv(sprintf('../%s/%s/%s/%s/common_age-sex_%s_%s_%s.tsv', status, cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
			pc10 = read.csv(sprintf('../%s/%s/%s/%s/common_pc10_%s_%s_%s.tsv', status, cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)

			k_filename_prefix = sprintf('common_kinship_%s_%s_%s', table_i, var_i, pheno_i)
			k = load_kinship(load_dir, k_filename_prefix)
			
			y=as.matrix(y)
			n = dim(y)[1]

			if (n != 0) {
				intercept = matrix(1, nrow = length(y))
				cov_df = cbind(intercept, var, age_sex, pc10)
				# cov_df[, 'sex'] <- as.numeric(as.factor(cov_df[, 'sex']))
				# cov_df[, 'gengrp6'] <- as.numeric(as.factor(cov_df[, 'gengrp6']))
				# cov_df[, 'CENTER'] <- as.numeric(as.factor(cov_df[, 'CENTER']))
				# cov_df[, 'WEIGHT_FINAL_NORM_OVERALL'] <- as.numeric(as.factor(cov_df[, 'WEIGHT_FINAL_NORM_OVERALL']))
				intercept_var_age_sex_pc10 = as.matrix(cov_df)
				cov_col = ncol(intercept_var_age_sex_pc10)
			}

			if (n == 0 | n <= cov_col) {
				print(sprintf("%s, %s too few samples including sample size = 0", var_i, pheno_i))
				beta_1 = NaN
				se_1 = NaN
				p_z_1 = NaN
			} else {
				k = as.matrix(k)
				print(sprintf("%s, %s started", var_i, pheno_i))
				if (pheno_i %in% biclass_list) {
					association = logistic.mm.aireml(Y=y, X=intercept_var_age_sex_pc10, k, verbose = FALSE, get.P = FALSE)
				} else {
					y = apply(y, 2, invnorm)
					association = lmm.aireml(Y=y, X=intercept_var_age_sex_pc10, k, verbose = FALSE, get.P = FALSE)
				}

				print(sprintf("%s, %s finished", var_i, pheno_i))
				beta_1 = association$BLUP_beta[2]
				se_1 = sqrt(association$varbeta[2, 2])
				score_1 = beta_1 / se_1
				p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
			}
			
			if (grepl("cn-binary_cn", var_i)) {
				n = sum(var[, 1])
			}
			
			if (grepl("-", var_i) && !grepl("cn-binary", var_i)) {
				table_df[pheno_i, sprintf("N_%s", var_i)] = sprintf("%d", n)
				table_df[pheno_i, sprintf("p-val_%s", var_i)] = sprintf("%s", p_z_1)
			} else {
				table_df[pheno_i, sprintf("N_%s", var_i)] = sprintf("%d", n)
				table_df[pheno_i, sprintf("beta_%s", var_i)] = sprintf("%s", beta_1)
				table_df[pheno_i, sprintf("SE_%s", var_i)] = sprintf("%s", se_1)
				table_df[pheno_i, sprintf("p-val_%s", var_i)] = sprintf("%s", p_z_1)
			}
		}
	}
}

write.table(table_df, file = sprintf("../%s/%s/%s_%s.tsv", status, cohort, cohort, table_i), row.names=TRUE, col.names=NA, sep='\t', quote=FALSE)
