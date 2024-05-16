library(tidyr)
library(dplyr)
library(purrr)
options(dplyr.summarise.inform = FALSE)


phenotype_suffix = snakemake@params[['phenotype_suffix']]
linear_model_res_path = snakemake@params[['linear_model_res_path']]
this_quantile = as.numeric(snakemake@params[['quantile']])
phenotypes = snakemake@params[['phenotypes']]
direction = 'top'

print(direction)
print(this_quantile)


print(paste('Using results from', linear_model_res_path))

compute_enrichment_zscore = function(df, zscore){
  enrich_df = df %>% filter(zscore_y > zscore) %>%
    group_by(model_name, model, burden, phenotype) %>%
    summarize(size = n(),
              extremes = sum(extreme_estimate)) %>%
    ungroup() %>%
    mutate(min_z_score = zscore)
  return(enrich_df)
}
z_score_thresholds = seq(0, 3, length.out = 1000)

res_list = list()
for (this_phenotype in phenotypes){
  print(this_phenotype)
  # this_phenotype = 'Triglycerides'
  
  this_df = readRDS(file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_sub_', this_phenotype, '_', phenotype_suffix,'.Rds'))) 
  print(nrow(this_df)/nrow(this_df %>% distinct(model)))
  
  quantile_df = this_df %>% #filter(gene_list == 'deeprvat_discoveries') %>%
    select(-gene_list) %>%
    group_by(model_name, model, burden,  phenotype) %>% 
    summarise(quant_thres = quantile(estimate, this_quantile)) %>%
    ungroup()
  
  
  input_df = this_df %>% 
    # drop_na(vtype_repeat) %>%
    select(-gene_list) %>%
    group_by(model_name, model, burden, phenotype) %>%
    mutate(zscore_y = (Y - mean(Y))/sd(Y)) %>%
    ungroup() %>%
    merge(quantile_df, how = 'left') %>%
    mutate(extreme_estimate = estimate >= quant_thres)
  
  if (direction == 'both'){
    print('taking absolute z-score to check enrichment in top and bottom extremes')
    input_df = input_df %>% mutate(zscore_y = abs(zscore_y))
  }
  
  print(max(input_df$zscore_y))
  z_score_thresholds = seq(0, 5.2, length.out = 1500)
  
  
  enrichment_df = map_dfr(z_score_thresholds, ~compute_enrichment_zscore(input_df, .x)) %>%
    mutate(prs_quantile = this_quantile)
  
  res_list = append(res_list, list(enrichment_df))
  
}
print('saving data')

direction_suffix = ifelse(direction == 'both', '_both', '')
combined_enrichment = do.call(rbind, res_list)
saveRDS(combined_enrichment, file.path(linear_model_res_path, 'plotting_data', paste0('enrichment_prs_vs_zscore_', direction_suffix,
                                                                                   this_quantile, '_', phenotype_suffix, '.Rds')))

