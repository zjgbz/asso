library(gaston)

args = commandArgs(trailingOnly=TRUE)
cohort = args[1]
#cohort = "JHS"
#cohort = "freeze6a-AA-noJHS"
#cohort = "SOL"
table_i = "table2"
var_list = list("rs334_del0", "rs334_del1", "rs334_del2", "del_rs334")
#pheno_list = list("hemoglobin_mcnc_bld_1", "hematocrit_vfr_bld_1", "rbc_ncnc_bld_1","mcv_entvol_rbc_1",
#                  "mch_entmass_rbc_1", "mchc_mcnc_rbc_1", "rdw_ratio_rbc_1", "lnHBA1C", "lnHBA1C_noDIABETES",
#                  "EGFRCKDEPI", "anemia", "microcytosis", "CKD")

pheno_list = list("hemoglobin_mcnc_bld_1", "hematocrit_vfr_bld_1", "rbc_ncnc_bld_1",
                  "mcv_entvol_rbc_1", "mch_entmass_rbc_1", "mchc_mcnc_rbc_1", "rdw_ratio_rbc_1",
                  "neutrophil_ncnc_bld_1", "lymphocyte_ncnc_bld_1", "basophil_ncnc_bld_1",
                  "eosinophil_ncnc_bld_1", "monocyte_ncnc_bld_1", "wbc_ncnc_bld_1",
                  "pmv_entvol_bld_1", "platelet_ncnc_bld_1", "EGFRCKDEPI", "CKD", "microcytosis", "anemia")

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
			if (var_i == "del_rs334") {
				var_inter = as.matrix(var$DEL * var$geno)
				intercept_var_age_sex_pc10 = as.matrix(cbind(intercept, var, var_inter, age_sex, pc10))
			} else {
				intercept_var_age_sex_pc10 = as.matrix(cbind(intercept, var, age_sex, pc10))
			}
			cov_col = ncol(intercept_var_age_sex_pc10)
		}

		if (n == 0 | n <= cov_col) {
			print(sprintf("%s, %s too few samples including sample size = 0", var_i, pheno_i))
			if (var_i == "del_rs334"){
				if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD"){
					p_z_1 = NaN
				} else {
					p_z_1 = NaN
				}
			
			} else {
				if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD"){
					or_1 = NaN
					ci_1_l = NaN
					ci_1_h = NaN
					p_z_1 = NaN
				} else {
					beta_1 = NaN
					se_1 = NaN
					p_z_1 = NaN
				}
			}
		} else {
			k = as.matrix(k)
			print(sprintf("%s, %s started", var_i, pheno_i))
			if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD") {
				association = logistic.mm.aireml(Y=y, X=intercept_var_age_sex_pc10, k, verbose = FALSE, get.P = FALSE)
			} else {
				association = lmm.aireml(Y=y, X=intercept_var_age_sex_pc10, k, verbose = FALSE, get.P = FALSE)
			}
			print(sprintf("%s, %s finished", var_i, pheno_i))

			if (var_i == "del_rs334"){
				if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD") {
					beta_1 = association$BLUP_beta[4]
					se_1 = sqrt(association$varbeta[4, 4])
					score_1 = beta_1 / se_1
					or_1 = exp(beta_1)
					p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
					ci_1_l = exp(beta_1 - 1.96 * se_1)
					ci_1_h = exp(beta_1 + 1.96 * se_1)
				} else {
					beta_1 = association$BLUP_beta[4]
					se_1 = sqrt(association$varbeta[4, 4])
					score_1 = beta_1 / se_1
					p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
				}
			} else {
				if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD") {
					beta_1 = association$BLUP_beta[2]
					se_1 = sqrt(association$varbeta[2, 2])
					score_1 = beta_1 / se_1
					or_1 = exp(beta_1)
					ci_1_l = exp(beta_1 - 1.96 * se_1)
					ci_1_h = exp(beta_1 + 1.96 * se_1)
					p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
				} else {
					beta_1 = association$BLUP_beta[2]
					se_1 = sqrt(association$varbeta[2, 2])
					score_1 = beta_1 / se_1
					p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
				}
			}
		}
		
		#p_t = 2 * pt(abs(score), df = n - 1, lower.tail = FALSE)
		
		sink(sprintf("../cohort/%s/%s.txt", cohort, table_i), append=TRUE)
		print(sprintf("%s, %s", var_i, pheno_i))
		print(n)
		if (var_i == "del_rs334"){
			if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD"){
				print(p_z_1)
			} else {
				print(p_z_1)
			}
		
		} else {
			if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD"){
				print(or_1)
				print(ci_1_l)
				print(ci_1_h)
				print(p_z_1)
			} else {
				print(beta_1)
				print(se_1)
				print(p_z_1)
			}
		}
		
		sink()
	}
}

