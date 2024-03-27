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
btypes = config[["btypes"]]


### run regression
combined_res = list()

# load y_data 
this_top_genes_to_keep = NULL
this_input_dir = file.path(regression_data_input_path)
log_info(this_input_dir)
y_data = LoadYData(this_input_dir)

# get phenotype name and corresponding covariate cols
all_column_names = GetColumnNames(phenotype_suffix, phenotype, covariates_dict)
covariate_cols = all_column_names[["covariate_cols"]]
phenotype_col = all_column_names[["phenotype_col"]]
log_info(paste("covariate_cols from GetColNames:", paste(covariate_cols, collapse = ", ")))


rare_burden_res = PredictPhenoBaseline(y_data, covariate_cols, phenotype_col, this_input_dir, top_q, gene_lists = gene_lists, btypes = btypes, genes_to_keep = this_top_genes_to_keep, use_top_q = use_top_q)

covarites_only_res = PredictPhenoCovariates(y_data, covariate_cols, phenotype_col, this_input_dir, top_q, use_top_q = use_top_q)


log_info("Combining data")

for (key in names(rare_burden_res)) {
    combined = gdata::combine(rare_burden_res[[key]], covarites_only_res[[key]], names = c("rare_burdens", "covariates_baseline")) %>%
        rename(method = source) %>%
        mutate(phenotype_col = phenotype_col, top_quantile = top_q)
    combined = rbind(combined_res[[key]], combined)
    combined_res[key] = list(combined)

}



log_info("writing result to {model_out_file}")
saveRDS(combined_res, model_out_file)
log_info("saved")
