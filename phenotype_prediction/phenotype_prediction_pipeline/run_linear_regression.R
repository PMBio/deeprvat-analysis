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
model_out_file = snakemake@params[["model_out_file"]]



regression_data_input_path = snakemake@params[["regression_data_input_path"]]
gene_lists = config[["gene_lists"]]
code_dir = config[["code_dir"]]

btypes = config[["btypes"]]

source(paste0(code_dir, "/regression_model_utils.R"))
source(paste0(code_dir, "/linear_regression_model_utils.R"))



################# Run pipeline from here ################################## ################# ################# #################

### run regression

combined_res = list()

this_top_genes_to_keep = NULL

# load y_data 
this_input_dir = file.path(regression_data_input_path)
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
log_info(paste("covariate_cols from GetColNames:", paste(covariate_cols, collapse = ", ")))

rare_burden_res = PredictPhenoBaselineLinear(y_data, covariate_cols, phenotype_col, this_input_dir, top_q, gene_lists = gene_lists, btypes = btypes, genes_to_keep = this_top_genes_to_keep)

covarites_only_res = PredictPhenoCovariatesLinear(y_data, covariate_cols, phenotype_col, this_input_dir, top_q)


for (key in names(rare_burden_res)) {
    combined = gdata::combine(rare_burden_res[[key]], covarites_only_res[[key]], names = c("rare_burdens", "covariates")) %>%
        rename(method = source) %>%
        mutate(phenotype_col = phenotype_col, top_quantile = top_q)
    combined = rbind(combined_res[[key]], combined)
    combined_res[key] = list(combined)

}

log_info("writing result to {model_out_file}")

saveRDS(combined_res, model_out_file)

