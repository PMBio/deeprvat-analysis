library(dplyr)
library(logger)
library(purrr)
library(tidyr)
library(stringr)

## params to specify #####

#######

phenotypes = snakemake@params[['phenotypes']]
linear_model_res_path = snakemake@params[['input_res_path']]
fdr = snakemake@params[['fdr']]
phenotype_suffix = snakemake@params[['phenotype_suffix']]
out_dir = snakemake@params[['out_dir']]
#######


phenotypes <- phenotypes[file.exists(paste0(linear_model_res_path, '/', phenotypes, '/', phenotype_suffix, '_fdr-', fdr, '.Rds'))]

p_list = list()
plot_df_list = list()

for (this_phenotype in phenotypes){
  print(this_phenotype)

  model_res_path = paste0(linear_model_res_path, '/', this_phenotype, '/', phenotype_suffix, '_fdr-', fdr, '.Rds')

  print(paste('Using model res path', model_res_path, sep = ' : '))
  combined_res = readRDS(model_res_path)
    
  res = combined_res$res
  
  res = res %>% filter(split == 'test') %>% select(-top_quantile, -estimate_bin, -truth, -split)
  
  print('Indexing samples according to order in data frame')
  IndexDF <- function(this_deeprvat_name){
    this_res = res %>% filter(model_name == this_deeprvat_name)
    this_res$index <- 1:nrow(this_res)
    return(this_res)
  }
  res = map_dfr(unique(res$model_name), IndexDF)
  print('Sanity check if rows with the same index have the same Y value')
  print(length(unique(res %>% filter(index == 1) %>% pull(Y))) == 1)
  
  print('computing deeprvat average predictions')
  deeprvat_avg_df = res %>% filter(model_name %in% grep('^deeprvat', unique(res$model_name), value = TRUE))
  print('extracting gene list from model name')
  deeprvat_avg_df[['gene_list']] = sapply(strsplit(deeprvat_avg_df[['model_name']], "-", fixed=TRUE), tail, 1)
  colnames(deeprvat_avg_df)
  deeprvat_avg_df = deeprvat_avg_df %>%
    group_by_at(vars(-estimate, -model_name)) %>% summarise(estimate = mean(estimate)) %>%
    ungroup() %>%
    mutate(model_name = paste('deeprvat-avg', gene_list, sep = '-')) %>%
    select(-gene_list)
  print('Preparing data for scatter plot')
  
  res = rbind(res, deeprvat_avg_df) %>% 
    arrange(model_name, index)  
  this_plot_df = res %>% #filter((grepl(this_gene_list,model_name) | grepl("covariates",model_name))) %>%
    filter(fdr == fdr)  %>% 
    separate(col = model_name, into = c('burden', 'vtype_repeat', 'gene_list'), sep = '-', remove = FALSE) %>%
    mutate(model = paste(burden, vtype_repeat, sep = '-')) %>%
    mutate(model = str_replace(model, "covariates-NA", "covariates")) %>%
    # filter(model %in% scatter_plot_models_to_keep) %>%
    mutate(phenotype = this_phenotype) 
  
  unique(this_plot_df$model_name)
  
  # this_plot_df[['model']] = factor(this_plot_df[['model']], levels = scatter_plot_models_to_keep)
  
  plot_df_list[[this_phenotype]] = this_plot_df %>% mutate(phenotype = this_phenotype)
  
}

saveRDS(plot_df_list, file.path(out_dir, paste0('plot_df_list_', phenotype_suffix, '_', fdr,'.Rds')))
