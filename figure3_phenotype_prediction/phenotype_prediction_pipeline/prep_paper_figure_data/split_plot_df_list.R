library(tidyr)
library(dplyr)

phenotype_suffix = snakemake@params[['phenotype_suffix']]
linear_model_res_path = snakemake@params[['linear_model_res_path']]


all_combined_metrics = readRDS(file.path(linear_model_res_path, 'plotting_data', paste0('all_recomputed_metrics_test_', phenotype_suffix,'.Rds')))

baseline_btypes = as.character(all_combined_metrics %>% distinct(model_name) %>% pull(model_name))
baseline_btypes =  unique(sapply(strsplit(baseline_btypes, "-"), function(x) x[1]))
baseline_btypes = grep('covariates|deeprvat', baseline_btypes, value = TRUE, invert = TRUE)

method_gene_list_filter_2 = quo((`burden`  %in% baseline_btypes & gene_list == 'baseline_only') | 
                                  (`burden` == 'deeprvat' & gene_list == 'deeprvat_discoveries'))
method_gene_list_filter_2

print('loading plot_df_list')

plot_df_list = readRDS(file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_', phenotype_suffix,'.Rds')))

for (p in names(plot_df_list)){
# for (p in c('LDL_direct', 'Triglycerides')){
  print(p)
  saveRDS(plot_df_list[[p]], file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_', p, '_', phenotype_suffix,'.Rds')))
  this_df_sub = plot_df_list[[p]] %>%
    select(estimate, Y, model_name, model, burden, gene_list, index,  phenotype) %>%
    filter(!!method_gene_list_filter_2 | (burden == 'covariates')) 
  # log_info('number of rows {nrow(this_df_sub)}')
  saveRDS(this_df_sub, file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_sub_', p, '_', phenotype_suffix,'.Rds')))
}
print('touching file')
file.create(file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_sub_', phenotype_suffix,'.finished')))