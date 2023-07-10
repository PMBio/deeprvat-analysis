library(Rcpp)
library(stringr)
library(arrow)
library(yardstick)
library(ggplot2)
library(gdata)
library(dplyr)
library(tidyr)
library(logger)

log_info("Snakemake params: {snakemake@params}")
config = snakemake@params[["config"]]
log_info("Config: {config}")
# use_rank = snakemake@params[['use_rank']]
use_rank = FALSE

phenotype = snakemake@params[["phenotype"]]
phenotype_suffix = snakemake@params[["phenotype_suffix"]]
# cv_split = snakemake@params[['cv_split']]
model_out_file = snakemake@params[["model_out_file"]]



regression_data_input_path = snakemake@params[["regression_data_input_path"]]
gene_lists = config[["gene_lists"]]
code_dir = config[["code_dir"]]

vtypes = config[["vtypes"]]
n_deeprvat_repeats = config[["n_deeprvat_repeats"]]
n_cv_splits = config[["n_cv_splits"]]

source(paste0(code_dir, "/regression_model_utils.R"))
source(paste0(code_dir, "/linear_regression_model_utils.R"))



################# Run pipeline from here ################################## ################# ################# #################

### run regression

combined_res = list()
for (cv_split in seq(0, n_cv_splits - 1)) {
    # seq(0, n_cv_splits-1)
    log_info("CV_split {cv_split},  phenotype_suffix: {phenotype_suffix}")

    this_top_genes_to_keep = NULL

    # load y_data which is constant for a given cv split
    this_input_dir = file.path(regression_data_input_path, paste("cv_split", cv_split, sep = ""))
    y_data = LoadYData(this_input_dir)
    if (!(any(grepl("_irnt$", colnames(y_data$train))))) {
        log_info("adding irnt column")
        for (key in names(y_data)) {
            y_data[[key]][[paste(phenotype, "irnt", sep = "_")]] = irnt(y_data[[key]][[phenotype]])
        }
    }

    top_q = 0.01  # just a dummy for LoadDataFor Model to generat Y_bin, which is not used by the linear model
    # get phenotype name and corresponding covariate cols
    all_column_names = GetColumnNames(phenotype_suffix, phenotype, covariates_dict)

    covariate_cols = all_column_names[["covariate_cols"]]
    phenotype_col = all_column_names[["phenotype_col"]]
    log_info("covariate_cols from GetColNames: {covariate_cols}")

    deeprvat_res = PredictPhenoDeepRVATLinear(y_data, covariate_cols, phenotype_col, this_input_dir, top_q, gene_lists = gene_lists, n_deeprvat_repeats = n_deeprvat_repeats, genes_to_keep = this_top_genes_to_keep)

    deeprvat_res_avg_burden = PredictPhenoDeepRVATLinearAvgBurden(y_data, covariate_cols, phenotype_col, this_input_dir, top_q, gene_lists = gene_lists, n_deeprvat_repeats = n_deeprvat_repeats, genes_to_keep = this_top_genes_to_keep)

    baseline_res = PredictPhenoBaselineLinear(y_data, covariate_cols, phenotype_col, this_input_dir, top_q, gene_lists = gene_lists, vtypes = vtypes, genes_to_keep = this_top_genes_to_keep)

    covarites_only_res = PredictPhenoCovariatesLinear(y_data, covariate_cols, phenotype_col, this_input_dir, top_q)


    for (key in names(deeprvat_res)) {
        combined = gdata::combine(deeprvat_res[[key]], deeprvat_res_avg_burden[[key]], baseline_res[[key]], covarites_only_res[[key]], names = c("deeprvat", "deeprvat_avg_burden", "rare_baseline", "covariates")) %>%
            rename(method = source) %>%
            mutate(cv_split = cv_split, phenotype_col = phenotype_col, top_quantile = top_q)
        combined = rbind(combined_res[[key]], combined)
        combined_res[key] = list(combined)

    }
}


saveRDS(combined_res, model_out_file)

