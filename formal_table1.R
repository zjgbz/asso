library(gaston)

args = commandArgs(trailingOnly=TRUE)
cohort = args[1]
#cohort = "JHS"
#cohort = "freeze6a-AA-noJHS"
#cohort = "SOL"
table_i = "table1"
biclass_list = list("CKD", "microcytosis", "anemia")
var_list = list("cn-binary_cn0", "cn-binary_cn1", "cn-binary_cn3", "rs334", "rs33930165")
if (cohort == "JHS") {
	pheno_list = list("hemoglobin_mcnc_bld_1", "hematocrit_vfr_bld_1", "rbc_ncnc_bld_1","mcv_entvol_rbc_1",
    	              "mch_entmass_rbc_1", "mchc_mcnc_rbc_1", "rdw_ratio_rbc_1", "lnHBA1C", "lnHBA1C_noDIABETES",
    	              "EGFRCKDEPI", "anemia", "microcytosis", "CKD")
} else {
	pheno_list = list("hemoglobin_mcnc_bld_1", "hematocrit_vfr_bld_1", "rbc_ncnc_bld_1",
	                  "mcv_entvol_rbc_1", "mch_entmass_rbc_1", "mchc_mcnc_rbc_1", "rdw_ratio_rbc_1",
	                  "neutrophil_ncnc_bld_1", "lymphocyte_ncnc_bld_1", "basophil_ncnc_bld_1",
	                  "eosinophil_ncnc_bld_1", "monocyte_ncnc_bld_1", "wbc_ncnc_bld_1",
	                  "pmv_entvol_bld_1", "platelet_ncnc_bld_1", "EGFRCKDEPI", "CKD", "microcytosis", "anemia")
}

header_list = list()
for (var_i in var_list) {
	header_list = append(header_list, sprintf("%s_N", var_i))
	header_list = append(header_list, sprintf("%s_beta", var_i))
	header_list = append(header_list, sprintf("%s_SE", var_i))
	header_list = append(header_list, sprintf("%s_p-val", var_i))
}

row_num = length(pheno_list)
col_num = length(header_list)
table_df = data.frame(matrix(ncol = col_num, nrow = row_num))
colnames(table_df) = header_list
row.names(table_df) = pheno_list

for (var_i in var_list) {
	for (pheno_i in pheno_list) {
		y = read.csv(sprintf('../cohort/%s/%s/%s/common_pheno_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		var = read.csv(sprintf('../cohort/%s/%s/%s/common_var_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		age_sex = read.csv(sprintf('../cohort/%s/%s/%s/common_age-sex_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		pc10 = read.csv(sprintf('../cohort/%s/%s/%s/common_pc10_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		k = read.csv(sprintf('../cohort/%s/%s/%s/common_kinship_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		
		y=as.matrix(y)
		n = dim(y)[1]

		if (n != 0) {
			intercept = matrix(1, nrow = length(y))
			cov_df = cbind(intercept, var, age_sex, pc10)
			cov_df[, 'sex'] <- as.numeric(as.factor(cov_df[, 'sex']))
			cov_df[, 'gengrp6'] <- as.numeric(as.factor(cov_df[, 'gengrp6']))
#			print(cov_df[, 'gengrp6'])
			cov_df[, 'CENTER'] <- as.numeric(as.factor(cov_df[, 'CENTER']))
			cov_df[, 'WEIGHT_FINAL_NORM_OVERALL'] <- as.numeric(as.factor(cov_df[, 'WEIGHT_FINAL_NORM_OVERALL']))
#			print(sapply(as.data.frame(cov_df), class))
			intercept_var_age_sex_pc10 = as.matrix(cov_df)
			
			cov_col = ncol(intercept_var_age_sex_pc10)
		}

		if (n == 0 | n <= cov_col) {
			print(sprintf("%s, %s too few samples including sample size = 0", var_i, pheno_i))
			if (pheno_i %in% biclass_list) {
				or_1 = NaN
				ci_1_l = NaN
				ci_1_h = NaN
				p_z_1 = NaN
				beta_1 = NaN
				se_1 = NaN
				p_z_1 = NaN
			} else {
				beta_1 = NaN
				se_1 = NaN
				p_z_1 = NaN
			}
		} else {
			k = as.matrix(k)
			print(sprintf("%s, %s started", var_i, pheno_i))
			if (pheno_i %in% biclass_list) {
				association = logistic.mm.aireml(Y=y, X=intercept_var_age_sex_pc10, k, verbose = FALSE, get.P = FALSE)
			} else {
				association = lmm.aireml(Y=y, X=intercept_var_age_sex_pc10, k, verbose = FALSE, get.P = FALSE)
			}

			print(sprintf("%s, %s finished", var_i, pheno_i))
			if (pheno_i %in% biclass_list) {
				beta_1_bi = association$BLUP_beta[2]
				se_1_bi = sqrt(association$varbeta[2, 2])
				score_1 = beta_1_bi / se_1_bi
				or_1 = exp(beta_1_bi)
				ci_1_l = exp(beta_1_bi - 1.96 * se_1_bi)
				ci_1_h = exp(beta_1_bi + 1.96 * se_1_bi)
				p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)

				beta_1 = exp(beta_1_bi)
				se_1 = sprintf("%.4f,%.4f", ci_1_l, ci_1_h)
				p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
			} else {
				beta_1 = association$BLUP_beta[2]
				se_1 = sqrt(association$varbeta[2, 2])
				score_1 = beta_1 / se_1
				p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
			}
		}
		
		table_df[pheno_i, sprintf("%s_N", var_i)] = sprintf("%d", n)
		table_df[pheno_i, sprintf("%s_beta", var_i)] = sprintf("%s", beta_1)
		table_df[pheno_i, sprintf("%s_SE", var_i)] = sprintf("%s", se_1)
		table_df[pheno_i, sprintf("%s_p-val", var_i)] = sprintf("%s", p_z_1)
	}
}

write.table(table_df, file = sprintf("../cohort/%s/%s_%s.tsv", cohort, cohort, table_i), row.names=TRUE, col.names=TRUE, sep='\t')
