library(dplyr)
library(logger)
library(purrr)
library(stringr)
library(yardstick)

## params to specify #####


phenotypes = snakemake@params[["phenotypes"]]
logistic_model_res_path = snakemake@params[["input_res_path"]]
phenotype_suffix = snakemake@params[["phenotype_suffix"]]
top_bottom_q_vals = snakemake@params[["top_bottom_q_vals"]]
top_quantiles = snakemake@params[["top_quantiles"]]
out_dir = snakemake@params[["out_dir"]]
code_dir = snakemake@params[["code_dir"]]
log_info("Analysing {length(phenotypes)} phenotypes")



##########################################
##########################################


GetRankedMetrics <- function(combined_res_ranked){
  
  recomputeMetrics <- function(this_res, group_info){
    top_bottom_q = group_info$extreme_direction
    decision_threshold = 0.5
    estimate = this_res[['estimate']]
    Y = this_res[['Y']]#unique(this_res$top_q)
    
    if (top_bottom_q == 'topq'){
      thres = quantile(Y, 1-unique(group_info$top_quantile))
      truth = as.integer(Y > thres)
    }else{
      thres = quantile(Y, unique(group_info$top_quantile))
      truth = as.integer(Y < thres)
    }
    fitted_results_bin <- ifelse(estimate > decision_threshold,1,0)
    this_res_new = tibble(estimate = estimate,
                          estimate_bin = as.factor(fitted_results_bin), 
                          truth = as.factor(truth),
                          Y = Y)
    
    class_and_probs_metrics <- metric_set(roc_auc, pr_auc, accuracy, f_meas)
    this_metrics = this_res_new %>% class_and_probs_metrics(truth, estimate,
                                                            estimate = estimate_bin, 
                                                            event_level = 'second')
    return(cbind(this_metrics, group_info))
  }
  
  combined_metrics_list = combined_res_ranked %>% 
    group_by(model_name, phenotype_col, top_quantile, extreme_direction) %>%
    group_map(recomputeMetrics)
  combined_metrics = do.call(rbind, combined_metrics_list) 
  return(combined_metrics)
}

combineModelResults <- function(phenotype, phenotype_suffix = NULL, top_q = NULL, model_res_dir,  top_bottom_q = "topq") {
    if (is.null(phenotype_suffix)) {
        this_res_files = list.files(paste(model_res_dir, phenotype, sep = "/"), full.names = TRUE)

    } else {
        this_res_files = list.files(model_res_dir, pattern = paste0(phenotype, '_', phenotype_suffix, "_", top_bottom_q, "-", top_q), full.names = TRUE)
    }
    this_res_files = this_res_files[grepl(top_bottom_q, this_res_files)]
    print(paste("this res files", this_res_files))

    # this_res_files = this_res_files[grep(fdr_rank_name, this_res_files)]
    # print(paste("this res files", this_res_files))
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
        mutate(extreme_direction = top_bottom_q)
    print(combined_res[, c("model_name", "method", "phenotype_col", "top_quantile")] %>%
        distinct())
    return(combined_res)
}


metrics_list = list()
for (phenotype in phenotypes) {
    for (top_q in top_quantiles) {
        for (top_bottom_q in top_bottom_q_vals) {
            log_info("getting combined res")
            log_info("{phenotype}, {top_q}, {top_bottom_q}")
            combined_res = combineModelResults(phenotype = phenotype, phenotype_suffix = phenotype_suffix, top_q = top_q, model_res_dir = logistic_model_res_path,  top_bottom_q = top_bottom_q)
            # res_list = append(res_list, list(combined_res %>% mutate(phenotype = phenotype)))

            metrics = GetRankedMetrics(combined_res)
            metrics_list = append(metrics_list, list(metrics %>%
                mutate(phenotype = phenotype)))

        }
    }
}

log_info("Combining and saving data")
combined_metrics = do.call(rbind, metrics_list)
log_info("Writing data to {out_dir}")

saveRDS(combined_metrics, file.path(out_dir, paste0("combined_metrics_", phenotype_suffix, ".Rds")))

log_info("Finished")




