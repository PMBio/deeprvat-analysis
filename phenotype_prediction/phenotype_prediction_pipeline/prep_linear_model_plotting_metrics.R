library(dplyr)
library(logger)
library(purrr)
library(yardstick)

## params to specify #####

phenotypes = snakemake@params[["phenotypes"]]
linear_model_res_path = snakemake@params[["input_res_path"]]
phenotype_suffix = snakemake@params[["phenotype_suffix"]]
out_dir = snakemake@params[["out_dir"]]

phenotypes <- phenotypes[file.exists(paste0(linear_model_res_path, "/", phenotypes,"_", phenotype_suffix, ".Rds"))]

recomputeContMetrics <- function(this_res, group_info) {
    this_metrics = metrics(this_res, Y, estimate)

    return(cbind(this_metrics, group_info))
}

RecomputeContMetrics <- function(phenotype) {
    print(phenotype)
    model_res_path = paste0(linear_model_res_path, "/", phenotype, "_", phenotype_suffix, ".Rds")
    # if (file.exists(model_res_path)){
    print(paste("Reading results from", model_res_path, sep = " : "))
    combined_res = readRDS(model_res_path)

    res = combined_res$res %>%
        select(-truth, -estimate_bin)

    metrics_list = res %>%
        group_by_at(vars(-estimate, -Y)) %>%
        group_map(recomputeContMetrics)
    combined_metrics = do.call(rbind, metrics_list)
    combined_metrics = combined_metrics %>%
        mutate_if(sapply(combined_metrics, is.character), as.factor) %>%
        mutate(phenotype = phenotype)
    return(combined_metrics)
}
all_combined_metrics = map_dfr(phenotypes, RecomputeContMetrics) %>%
    select(-.estimator)
saveRDS(all_combined_metrics, file.path(out_dir, paste0("all_recomputed_metrics_test_", phenotype_suffix, ".Rds")))




