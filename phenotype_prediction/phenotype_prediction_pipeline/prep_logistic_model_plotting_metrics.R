library(dplyr)
library(logger)
library(purrr)
library(stringr)

## params to specify #####


phenotypes = snakemake@params[["phenotypes"]]
logistic_model_res_path = snakemake@params[["input_res_path"]]
fdr = snakemake@params[["fdr"]]
phenotype_suffix = snakemake@params[["phenotype_suffix"]]
top_bottom_q_vals = snakemake@params[["top_bottom_q_vals"]]
top_quantiles = snakemake@params[["top_quantiles"]]
out_dir = snakemake@params[["out_dir"]]
code_dir = snakemake@params[["code_dir"]]
######
phenotypes <- phenotypes[dir.exists(paste0(logistic_model_res_path, "/", phenotypes))]
log_info("Analysing {length(phenotypes)} phenotypes")



combineModelResults <- function(phenotype, phenotype_suffix = NULL, top_q = NULL, fdr = 0.05, model_res_dir, use_rank = TRUE, top_bottom_q = "topq") {
    fdr_rank_name = ifelse(use_rank, "rank", "fdr")
    if (is.null(phenotype_suffix)) {
        this_res_files = list.files(paste(model_res_dir, phenotype, sep = "/"), full.names = TRUE)

    } else {
        this_res_files = list.files(paste(model_res_dir, phenotype, sep = "/"), pattern = paste0(phenotype_suffix, "_", top_bottom_q, "-", top_q, "_", fdr_rank_name, "-", fdr), full.names = TRUE)
    }
    this_res_files = this_res_files[grepl(top_bottom_q, this_res_files)]
    print(paste("this res files", this_res_files))

    this_res_files = this_res_files[grep(fdr_rank_name, this_res_files)]
    print(paste("this res files", this_res_files))
    print(length(this_res_files))
    res_list = list()
    for (file in this_res_files) {
        print(file)
        this_res = readRDS(file)[["res"]]
        res_list = append(res_list, list(this_res))
    }
    print("combining data")
    combined_res = do.call(rbind, res_list)
    combined_res = combined_res %>%
        mutate(extreme_direction = top_bottom_q, `:=`(!!fdr_rank_name, fdr))
    print(combined_res[, c(fdr_rank_name, "model_name", "method", "phenotype_col", "top_quantile")] %>%
        distinct())
    return(combined_res)
}


metrics_list = list()
for (phenotype in phenotypes) {
    for (top_q in top_quantiles) {
        for (top_bottom_q in top_bottom_q_vals) {
            log_info("getting combined res")
            log_info("{phenotype}, {top_q}, {top_bottom_q}")
            combined_res = combineModelResults(phenotype = phenotype, phenotype_suffix = phenotype_suffix, top_q = top_q, model_res_dir = logistic_model_res_path, fdr = fdr, use_rank = FALSE, top_bottom_q = top_bottom_q)
            # res_list = append(res_list, list(combined_res %>% mutate(phenotype = phenotype)))

            metrics = GetRankedMetrics(combined_res, fdr_rank_name = "fdr")
            metrics_list = append(metrics_list, list(metrics %>%
                mutate(phenotype = phenotype)))

        }
    }
}

log_info("Combining and saving data")
combined_metrics = do.call(rbind, metrics_list)
# combined_res = do.call(rbind, res_list)
log_info("Writing data to {out_dir}")

# saveRDS(combined_res, file.path(out_dir, paste0('combined_res_', phenotype_suffix, '_', fdr, '.Rds')))
saveRDS(combined_metrics, file.path(out_dir, paste0("combined_metrics_", phenotype_suffix, "_", fdr, ".Rds")))

log_info("Finished")




