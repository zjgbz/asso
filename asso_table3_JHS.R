library(gaston)

cohort = "JHS"
table_i = "table3"
var_list = list("rs11248850", "rs11248850_del0", "rs11248850_del12", "del_rs11248850")
pheno_list = list("hemoglobin_mcnc_bld_1", "hematocrit_vfr_bld_1", "rbc_ncnc_bld_1",
	"mcv_entvol_rbc_1","mch_entmass_rbc_1", "mchc_mcnc_rbc_1", "rdw_ratio_rbc_1")

for (var_i in var_list) {
	for (pheno_i in pheno_list) {
		y = read.csv(sprintf('../cohort/%s/%s/%s/common_pheno_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		var = read.csv(sprintf('../cohort/%s/%s/%s/common_var_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		age_sex = read.csv(sprintf('../cohort/%s/%s/%s/common_age-sex_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		pc10 = read.csv(sprintf('../cohort/%s/%s/%s/common_pc10_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
		k = read.csv(sprintf('../cohort/%s/%s/%s/common_kinship_%s_%s_%s.tsv', cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)

		y=as.matrix(y)
		intercept = matrix(1, nrow = length(y))

		if (var_i == "del_rs11248850") {
			var_inter = as.matrix(var$DEL * var$geno)
			intercept_var_age_sex_pc10 = as.matrix(cbind(intercept, var, var_inter, age_sex, pc10))
		} else {
			intercept_var_age_sex_pc10 = as.matrix(cbind(intercept, var, age_sex, pc10))
		}
		
		k = as.matrix(k)
		n = dim(y)[1]
		
		print(sprintf("%s, %s started", var_i, pheno_i))
		association = lmm.aireml(Y=y, X=intercept_var_age_sex_pc10, k, verbose = FALSE, get.P = FALSE)
		print(sprintf("%s, %s finished", var_i, pheno_i))
		
		if (var_i == "del_rs11248850"){
			beta_1 = association$BLUP_beta[4]
			se_1 = sqrt(association$varbeta[4, 4])
			score_1 = beta_1 / se_1
			p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)		
		} else {
			beta_1 = association$BLUP_beta[2]
			se_1 = sqrt(association$varbeta[2, 2])
			score_1 = beta_1 / se_1
			p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
		}
		
		#p_t = 2 * pt(abs(score), df = n - 1, lower.tail = FALSE)
		
		sink(sprintf("../cohort/%s/%s.txt", cohort, table_i), append=TRUE)
		print(sprintf("%s, %s", var_i, pheno_i))
		print(n)
		if (var_i == "del_rs11248850"){
			print(p_z_1)
		} else {
			print(beta_1)
			print(se_1)
			print(p_z_1)
		}
		
		sink()
	}
}

