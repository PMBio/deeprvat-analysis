library(tibble)
library(stringr)
library(arrow)
library(yardstick)
library(ggplot2)
library(gdata)
library(dplyr)
library(logger)

config = snakemake@params[["config"]]
# use_rank = snakemake@params[['use_rank']]

phenotype = snakemake@params[["phenotype"]]
top_q = as.numeric(snakemake@params[["top_q"]])


phenotype_suffix = snakemake@params[["phenotype_suffix"]]
# cv_split = snakemake@params[['cv_split']]
log_info("model config {config}")
model_out_file = snakemake@params[["model_out_file"]]


use_top_q_dict = c(TRUE, FALSE)
names(use_top_q_dict) = c("topq", "bottomq")
top_bottom_q = snakemake@params[["top_bottom_q"]]
use_top_q = ifelse(use_top_q_dict[[top_bottom_q]], TRUE, FALSE)

regression_data_input_path = snakemake@params[["regression_data_input_path"]]
gene_lists = config[["gene_lists"]]
code_dir = config[["code_dir"]]

source(file.path(code_dir, "regression_model_utils.R"))

# other parameters
vtypes = config[["vtypes"]]
n_deeprvat_repeats = config[["n_deeprvat_repeats"]]
n_cv_splits = config[["n_cv_splits"]]


### run regression
combined_res = list()
for (cv_split in seq(0, n_cv_splits - 1)) {
    log_info("cv_split: {cv_split}  top_q: {top_q}, phenotype_suffix: {phenotype_suffix}")

    # load y_data which is constant for a given cv split
    this_top_genes_to_keep = NULL
    this_input_dir = file.path(regression_data_input_path, paste("cv_split", cv_split, sep = ""))
    log_info(this_input_dir)
    y_data = LoadYData(this_input_dir)

    # get phenotype name and corresponding covariate cols
    all_column_names = GetColumnNames(phenotype_suffix, phenotype, covariates_dict)
    covariate_cols = all_column_names[["covariate_cols"]]
    phenotype_col = all_column_names[["phenotype_col"]]
    log_info("covariate_cols from GetColNames {covariate_cols}")

    deeprvat_res = PredictPhenoDeepRVAT(y_data, covariate_cols, phenotype_col, this_input_dir, top_q, gene_lists = gene_lists, n_deeprvat_repeats = n_deeprvat_repeats, genes_to_keep = this_top_genes_to_keep, use_top_q = use_top_q)

    baseline_res = PredictPhenoBaseline(y_data, covariate_cols, phenotype_col, this_input_dir, top_q, gene_lists = gene_lists, vtypes = vtypes, genes_to_keep = this_top_genes_to_keep, use_top_q = use_top_q)

    covarites_only_res = PredictPhenoCovariates(y_data, covariate_cols, phenotype_col, this_input_dir, top_q, use_top_q = use_top_q)


    log_info("Combining data")

    for (key in names(deeprvat_res)) {
        combined = gdata::combine(deeprvat_res[[key]], baseline_res[[key]], covarites_only_res[[key]], names = c("deeprvat", "rare_baseline", "covariates_baseline")) %>%
            rename(method = source) %>%
            mutate(cv_split = cv_split, phenotype_col = phenotype_col, top_quantile = top_q)
        combined = rbind(combined_res[[key]], combined)
        combined_res[key] = list(combined)

    }
}


log_info("writing result to {model_out_file}")
saveRDS(combined_res, model_out_file)
log_info("saved")
