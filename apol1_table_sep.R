library(gaston)
library(feather)

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

gaston_logi_result <- function(gaston_model, var_col_number) {
	beta_1_bi = gaston_model$BLUP_beta[var_col_number]
	se_1_bi = sqrt(gaston_model$varbeta[var_col_number, var_col_number])
	score_1 = beta_1_bi / se_1_bi
	or_1 = exp(beta_1_bi)
	ci_1_l = exp(beta_1_bi - 1.96 * se_1_bi)
	ci_1_h = exp(beta_1_bi + 1.96 * se_1_bi)
	p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
	
	beta_1 = exp(beta_1_bi)
	se_1 = sprintf("%.4f,%.4f", ci_1_l, ci_1_h)
	p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
	gaston_result = list("beta" = beta_1, "se" = se_1, "p_val" = p_z_1)
	return (gaston_result)
}

gaston_lm_result <- function(gaston_model, var_col_number) {
	beta_1 = gaston_model$BLUP_beta[var_col_number]
	se_1 = sqrt(gaston_model$varbeta[var_col_number, var_col_number])
	score_1 = beta_1 / se_1
	p_z_1 = 2 * pnorm(abs(score_1), lower.tail = FALSE)
	gaston_result = list("beta" = beta_1, "se" = se_1, "p_val" = p_z_1)
	return (gaston_result)
}

args = commandArgs(trailingOnly=TRUE)
status_name = args[1]
cohort = args[2]
var_list_idx = strtoi(args[3])
pheno = args[4]
expand_interaction = args[5]

table_i_vec = c("rs2302524", "rs2302524", "rs2633317", "rs2633317", "rs4251805", "rs4251805",
	            "rs4760", "rs4760", "rs73935023", "rs73935023")
biclass_list = list("CKD", "microcytosis", "anemia")

var_list_list = list(list("rs2302524"), list("apol1-rs2302524"), list("rs2633317"), list("apol1-rs2633317"),
					 list("rs4251805"), list("apol1-rs4251805"), list("rs4760"), list("apol1-rs4760"),
					 list("rs73935023"), list("apol1-rs73935023"))
var_list = var_list_list[[var_list_idx]]
table_i = table_i_vec[var_list_idx]

# pheno_list_list = list(list("EGFRCKDEPI"), list("CKD"))
# pheno_list = pheno_list_list[[strtoi(args[3])]]
pheno_list = list(pheno)

header_list = list()
for (var_i in var_list) {
	header_list = append(header_list, sprintf("N_%s", var_i))
	header_list = append(header_list, sprintf("beta_%s", var_i))
	header_list = append(header_list, sprintf("SE_%s", var_i))
	header_list = append(header_list, sprintf("p-val_%s", var_i))

	if (expand_interaction == "expand_int" && (grepl("-", var_i) && !grepl("cn-binary", var_i))) {
		var_i_split = strsplit(var_i, "-")
		inter1 = var_i_split[[1]][1]
		inter2 = var_i_split[[1]][2]

		header_list = append(header_list, sprintf("beta_%s", inter1))
		header_list = append(header_list, sprintf("SE_%s", inter1))
		header_list = append(header_list, sprintf("p-val_%s", inter1))

		header_list = append(header_list, sprintf("beta_%s", inter2))
		header_list = append(header_list, sprintf("SE_%s", inter2))
		header_list = append(header_list, sprintf("p-val_%s", inter2))
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
		load_dir = file.path('..', status_name, cohort, table_i, var_i)
		y_dir_filename = sprintf('../%s/%s/%s/%s/common_pheno_%s_%s_%s.tsv', status_name, cohort, table_i, var_i, table_i, var_i, pheno_i)
		if (!(file.exists(y_dir_filename))) {
			next
		} else if (file.exists(y_dir_filename)) {
			y = read.csv(sprintf('../%s/%s/%s/%s/common_pheno_%s_%s_%s.tsv', status_name, cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
			var = read.csv(sprintf('../%s/%s/%s/%s/common_var_%s_%s_%s.tsv', status_name, cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
			age_sex = read.csv(sprintf('../%s/%s/%s/%s/common_age-sex_%s_%s_%s.tsv', status_name, cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)
			pc10 = read.csv(sprintf('../%s/%s/%s/%s/common_pc10_%s_%s_%s.tsv', status_name, cohort, table_i, var_i, table_i, var_i, pheno_i), sep='\t', header=TRUE, row.names=1)

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
					y = apply(y, 2, invnorm)
					association = lmm.aireml(Y=y, X=intercept_var_age_sex_pc10, k, verbose = FALSE, get.P = FALSE)
				}

				print(sprintf("%s, %s finished", var_i, pheno_i))
				if (pheno_i %in% biclass_list) {
					logi_result_1 = gaston_logi_result(association, 2)
					beta_1 = logi_result_1$beta
					se_1 = logi_result_1$se
					p_z_1 = logi_result_1$p_val
					if (expand_interaction == "expand_int" && (grepl("-", var_i) && !grepl("cn-binary", var_i))) {
						logi_result_2 = gaston_logi_result(association, 3)
						logi_result_3 = gaston_logi_result(association, 4)
						
						beta_2 = logi_result_2$beta
						se_2 = logi_result_2$se
						p_z_2 = logi_result_2$p_val

						beta_3 = logi_result_3$beta
						se_3 = logi_result_3$se
						p_z_3 = logi_result_3$p_val
					}
				} else {
					lm_result_1 = gaston_lm_result(association, 2)
					beta_1 = lm_result_1$beta
					se_1 = lm_result_1$se
					p_z_1 = lm_result_1$p_val
					if (expand_interaction == "expand_int" && (grepl("-", var_i) && !grepl("cn-binary", var_i))) {
						lm_result_2 = gaston_lm_result(association, 3)
						lm_result_3 = gaston_lm_result(association, 4)
						
						beta_2 = lm_result_2$beta
						se_2 = lm_result_2$se
						p_z_2 = lm_result_2$p_val

						beta_3 = lm_result_3$beta
						se_3 = lm_result_3$se
						p_z_3 = lm_result_3$p_val
					}
				}
			}
			
			table_df[pheno_i, sprintf("N_%s", var_i)] = sprintf("%d", n)
			table_df[pheno_i, sprintf("beta_%s", var_i)] = sprintf("%s", beta_1)
			table_df[pheno_i, sprintf("SE_%s", var_i)] = sprintf("%s", se_1)
			table_df[pheno_i, sprintf("p-val_%s", var_i)] = sprintf("%s", p_z_1)
			if (expand_interaction == "expand_int" && (grepl("-", var_i) && !grepl("cn-binary", var_i))) {
				var_i_split = strsplit(var_i, "-")
				inter1 = var_i_split[[1]][1]
				inter2 = var_i_split[[1]][2]
				
				table_df[pheno_i, sprintf("beta_%s", inter1)] = sprintf("%s", beta_2)
				table_df[pheno_i, sprintf("SE_%s", inter1)] = sprintf("%s", se_2)
				table_df[pheno_i, sprintf("p-val_%s", inter1)] = sprintf("%s", p_z_2)

				table_df[pheno_i, sprintf("beta_%s", inter2)] = sprintf("%s", beta_3)
				table_df[pheno_i, sprintf("SE_%s", inter2)] = sprintf("%s", se_3)
				table_df[pheno_i, sprintf("p-val_%s", inter2)] = sprintf("%s", p_z_3)
			}
		}
	}
}

write.table(table_df, file = sprintf("../%s/%s/%s_%s_%s.tsv", status_name, cohort, table_i, args[2], args[3]), row.names=TRUE, col.names=TRUE, sep='\t')
