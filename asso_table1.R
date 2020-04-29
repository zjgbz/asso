library(gaston)

args = commandArgs(trailingOnly=TRUE)
cohort = args[1]
#cohort = "JHS"
#cohort = "freeze6a-AA-noJHS"
#cohort = "SOL"
table_i = "table1"
var_list = list("cn-binary", "rs334", "rs33930165")
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
#		print("y start.")
		y = read.csv(sprintf('../cohort/%s/%s/%s/common_pheno_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
#		print("var start.")
		var = read.csv(sprintf('../cohort/%s/%s/%s/common_var_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
#		print("age_sex start.")
#		print(dim(var))
		age_sex = read.csv(sprintf('../cohort/%s/%s/%s/common_age-sex_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
#		print("pc10 start.")
#		print(dim(age_sex))
		pc10 = read.csv(sprintf('../cohort/%s/%s/%s/common_pc10_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
#		print("k start.")
#		print(dim(pc10))
		k = read.csv(sprintf('../cohort/%s/%s/%s/common_kinship_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		
#		print("y start.")
		y=as.matrix(y)
#		print(y)
		n = dim(y)[1]
#		print(n)

		if (n != 0) {
			intercept = matrix(1, nrow = length(y))
			intercept_var_age_sex_pc10 = as.matrix(cbind(intercept, var, age_sex, pc10))
			cov_col = ncol(intercept_var_age_sex_pc10)
		}

		if (n == 0 | n <= cov_col) {
			print(sprintf("%s, %s too few samples including sample size = 0", var_i, pheno_i))
			if (var_i == "cn-binary"){
				if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD"){
					or_1 = NaN
					ci_1_l = NaN
					ci_1_h = NaN
					p_z_1 = NaN
					or_2 = NaN
					ci_2_l = NaN
					ci_2_h = NaN
					p_z_2 = NaN
				} else {
					beta_1 = NaN
					se_1 = NaN
					p_z_1 = NaN
					beta_2 = NaN
					se_2 = NaN
					p_z_2 = NaN
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
			if (var_i == "cn-binary"){
				if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD") {
					beta_1 = association$BLUP_beta[3]
					se_1 = sqrt(association$varbeta[3, 3])
					score_1 = beta_1 / se_1
					beta_2 = association$BLUP_beta[2]
					se_2 = sqrt(association$varbeta[2, 2])
					score_2 = beta_2 / se_2
					or_1 = exp(beta_1)
					or_2 = exp(beta_2)
					p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
					p_z_2 = 2 * pnorm(abs(score_2), lower.tail = FALSE)
					ci_1_l = exp(beta_1 - 1.96 * se_1)
					ci_1_h = exp(beta_1 + 1.96 * se_1)
					ci_2_l = exp(beta_2 - 1.96 * se_2)
					ci_2_h = exp(beta_2 + 1.96 * se_2)
				} else {
					beta_1 = association$BLUP_beta[3]
					se_1 = sqrt(association$varbeta[3, 3])
					score_1 = beta_1 / se_1
					beta_2 = association$BLUP_beta[2]
					se_2 = sqrt(association$varbeta[2, 2])
					score_2 = beta_2 / se_2
					p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
					p_z_2 = 2 * pnorm(abs(score_2), lower.tail = FALSE)
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
		if (var_i == "cn-binary"){
			if (pheno_i == "anemia" || pheno_i == "microcytosis" || pheno_i == "CKD"){
				print(or_1)
				print(ci_1_l)
				print(ci_1_h)
				print(p_z_1)
				print(or_2)
				print(ci_2_l)
				print(ci_2_h)
				print(p_z_2)
			} else {
				print(beta_1)
				print(se_1)
				print(p_z_1)
				print(beta_2)
				print(se_2)
				print(p_z_2)
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

