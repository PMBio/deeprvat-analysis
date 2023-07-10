library(tidyr)
library(dplyr)
library(purrr)



exp_name = snakemake@params[['exp_name']]
fdr = snakemake@params[['fdr']]
phenotype_suffix = snakemake@params[['phenotype_suffix']]
linear_model_res_path = snakemake@params[['linear_model_res_path']]
outlier_mode = snakemake@params[['outlier_mode']]
extreme_quantile = as.numeric(snakemake@params[['extreme_quantile']])
phenotypes = snakemake@params[['phenotypes']]

print(extreme_quantile)
print(outlier_mode)


print(paste('Using results from', linear_model_res_path))
p_list_delta = list()

ranked_df_list = list()
all_ranks = seq(1, 5000, 10)


for (this_phenotype in phenotypes){
  print(this_phenotype)
  this_df = readRDS(file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_sub_', this_phenotype, '_', phenotype_suffix, '_', fdr,'.Rds'))) %>%
    select(-gene_list)
  unique(this_df$gene_list)
  sd_prs = sd(this_df %>% filter(model == 'covariates-PRS') %>% pull(estimate))
  
  #dummy rows for covariates for all gene lists
  covariates_res = this_df %>% filter(grepl('covariates', model_name))
  for (this_gene_list in na.omit(unique(this_df$gene_list))){
    this_df = rbind(this_df, covariates_res %>% mutate(gene_list = this_gene_list))
  }
  plot_df_delta = this_df  %>%
    select(model, estimate, Y, index) %>%
    pivot_wider(names_from = model, values_from = estimate)
  
  plot_df_delta = plot_df_delta %>%
    pivot_longer(grep('baseline|deeprvat', names(plot_df_delta), value = TRUE),
                 names_to = "Rare estimator", values_to = 'rare_prediction') %>%
    mutate(d_rare_prs = abs(`covariates-PRS` - rare_prediction))
  
  # define extreme individuals base on their phenotype
  all_y = plot_df_delta %>% select(Y, index) %>% distinct() %>% pull(Y)
  print(extreme_quantile)
  y_quantiles = quantile(all_y, probs = c(extreme_quantile, 1-extreme_quantile))
  if (outlier_mode == 'both'){
    half_extreme_quantile = extreme_quantile * 0.5
    y_quantiles_half = quantile(all_y, probs = c(half_extreme_quantile, 1-half_extreme_quantile))
    print(paste0('using top and bottom quantile'))
    plot_df_delta = plot_df_delta %>%
      mutate(is_outlier_measurment = ((Y < y_quantiles_half[[1]]) | (Y > y_quantiles_half[[2]])))
  }else if (outlier_mode == 'topq'){
    print(paste0('using top quantile >', y_quantiles[[2]]))
    plot_df_delta = plot_df_delta %>%
      mutate(is_outlier_measurment = (Y > y_quantiles[[2]]))
  }else if (outlier_mode == 'bottomq'){
    print(paste0('using bottom quantile <', y_quantiles[[1]]))
    plot_df_delta = plot_df_delta %>%
      mutate(is_outlier_measurment = (Y < y_quantiles[[1]]))
  }else{
    print('invalid outlier mode')
  }
  
  plot_df_delta = plot_df_delta %>%
    mutate(d_rare_y = abs(Y - rare_prediction), d_prs_y = abs(Y - `covariates-PRS`)) %>%
    mutate(is_improved_rare = d_rare_y < d_prs_y)
  
  print('computing metric for each rank')
  compute_rank_matrics = function(df, max_rank){
    ranked_df = plot_df_delta  %>%
      arrange(desc(d_rare_prs)) %>%
      group_by(`Rare estimator`) %>%
      mutate(rank=row_number())
    ranked_df = ranked_df %>% filter(rank < max_rank) %>%
      group_by(`Rare estimator`) %>%
      summarise(sum = sum(is_outlier_measurment)) %>% ungroup() %>%
      mutate(max_rank = max_rank)
    return(ranked_df)
  }
  ranked_df = map_dfr(all_ranks, ~compute_rank_matrics(plot_df_delta, .x)) %>%
    mutate(phenotype = this_phenotype) %>%
    mutate(extreme_quantile = extreme_quantile) %>%
    mutate(outlier_mode = outlier_mode)
  ranked_df_list = append(ranked_df_list, list(ranked_df))
  
}


print('Combining and saving data')

ranked_df_combined = do.call(rbind, ranked_df_list)
saveRDS(ranked_df_combined, file.path(linear_model_res_path, 'plotting_data', 'replication_in_extreme', paste0('ranked_df_combined_', 
                                                    outlier_mode, '_', extreme_quantile, '_', phenotype_suffix, '_', fdr,'.Rds')))

