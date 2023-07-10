library(tidyr)
library(dplyr)

fdr = snakemake@params[['fdr']]
phenotype_suffix = snakemake@params[['phenotype_suffix']]
linear_model_res_path = snakemake@params[['linear_model_res_path']]


all_combined_metrics = readRDS(file.path(linear_model_res_path, 'plotting_data', paste0('all_recomputed_metrics_test_', phenotype_suffix, '_', fdr,'.Rds')))

baseline_btypes = unique(sapply(strsplit(grep('baseline-', all_combined_metrics %>% distinct(model_name) %>% pull(model_name), value = TRUE), "-"), function(x) x[2]))

baseline_methods =  paste('baseline', baseline_btypes, sep = '-')

method_gene_list_filter_2 = quo((`burden` == 'baseline' & `vtype_repeat` %in% baseline_btypes & gene_list == 'baseline_only') | 
                                  (`burden` == 'deeprvat' & `vtype_repeat` == 'avg' & gene_list == 'deeprvat_discoveries'))
method_gene_list_filter_2

print('loading plot_df_list')

plot_df_list = readRDS(file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_', phenotype_suffix, '_', fdr,'.Rds')))

for (p in names(plot_df_list)){
# for (p in c('LDL_direct', 'Triglycerides')){
  print(p)
  saveRDS(plot_df_list[[p]], file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_', p, '_', phenotype_suffix, '_', fdr,'.Rds')))
  this_df_sub = plot_df_list[[p]] %>%
    select(estimate, Y, model_name, model, burden, vtype_repeat, gene_list, cv_split, index, phenotype) %>%
    filter(!!method_gene_list_filter_2 | (burden == 'covariates')) %>%
    filter(!grepl('_', vtype_repeat)) 
  
  saveRDS(this_df_sub, file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_sub_', p, '_', phenotype_suffix, '_', fdr,'.Rds')))
}
print('touching file')
file.create(file.path(linear_model_res_path, 'plotting_data', paste0('plot_df_list_sub_', phenotype_suffix, '_', fdr,'.finished')))