library(dplyr)
library(logger)
library(purrr)
library(yardstick)

## params to specify #####


phenotypes = snakemake@params[["phenotypes"]]
linear_model_res_path = snakemake@params[["input_res_path"]]
fdr = snakemake@params[["fdr"]]
phenotype_suffix = snakemake@params[["phenotype_suffix"]]
out_dir = snakemake@params[["out_dir"]]
#######

phenotypes <- phenotypes[file.exists(paste0(linear_model_res_path, "/", phenotypes, "/", phenotype_suffix, "_fdr-", fdr, ".Rds"))]

recomputeContMetrics <- function(this_res, group_info) {
    this_metrics = metrics(this_res, Y, estimate)

    return(cbind(this_metrics, group_info))
}

RecomputeContMetrics <- function(phenotype) {
    print(phenotype)
    model_res_path = paste0(linear_model_res_path, "/", phenotype, "/", phenotype_suffix, "_fdr-", fdr, ".Rds")
    # if (file.exists(model_res_path)){
    print(paste("Reading results from", model_res_path, sep = " : "))
    combined_res = readRDS(model_res_path)

    res = combined_res$res %>%
        filter(split == "test") %>%
        select(-truth, -estimate_bin, -split)

    print("computing deeprvat average predictions")
    IndexDF <- function(this_deeprvat_name) {
        this_res = res %>%
            filter(model_name == this_deeprvat_name)
        this_res$index <- 1:nrow(this_res)
        return(this_res)
    }
    deeprvat_avg_df = map_dfr(grep("^deeprvat", unique(res$model_name), value = TRUE), IndexDF)
    deeprvat_avg_df[["gene_list"]] = sapply(strsplit(deeprvat_avg_df[["model_name"]], "-", fixed = TRUE), tail, 1)
    colnames(deeprvat_avg_df)
    deeprvat_avg_df = deeprvat_avg_df %>%
        group_by_at(vars(-estimate, -model_name)) %>%
        summarise(estimate = mean(estimate)) %>%
        ungroup() %>%
        mutate(model_name = paste("deeprvat-avg", gene_list, sep = "-")) %>%
        select(-index, -gene_list)

    print("computing metrics")

    res = rbind(res, deeprvat_avg_df)
    metrics_list = res %>%
        group_by_at(vars(-estimate, -Y, -cv_split)) %>%
        group_map(recomputeContMetrics)
    combined_metrics = do.call(rbind, metrics_list)
    combined_metrics = combined_metrics %>%
        mutate_if(sapply(combined_metrics, is.character), as.factor) %>%
        mutate(phenotype = phenotype)
    return(combined_metrics)
}
all_combined_metrics = map_dfr(phenotypes, RecomputeContMetrics) %>%
    select(-.estimator)
saveRDS(all_combined_metrics, file.path(out_dir, paste0("all_recomputed_metrics_test_", phenotype_suffix, "_", fdr, ".Rds")))




