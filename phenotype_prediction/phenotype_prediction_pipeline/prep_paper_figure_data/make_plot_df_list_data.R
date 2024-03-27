library(dplyr)
library(logger)
library(purrr)
library(tidyr)
library(stringr)

## params to specify #####

#######

phenotypes = snakemake@params[['phenotypes']]
linear_model_res_path = snakemake@params[['input_res_path']]
phenotype_suffix = snakemake@params[['phenotype_suffix']]
out_dir = snakemake@params[['out_dir']]
#######


phenotypes <- phenotypes[file.exists(paste0(linear_model_res_path, '/', phenotypes, '_', phenotype_suffix, '.Rds'))]

p_list = list()
plot_df_list = list()

for (this_phenotype in phenotypes){
  print(this_phenotype)

  model_res_path = paste0(linear_model_res_path, '/', this_phenotype, '_', phenotype_suffix, '.Rds')

  print(paste('Using model res path', model_res_path, sep = ' : '))
  combined_res = readRDS(model_res_path)
    
  res = combined_res$res
  
  res = res %>% select(-top_quantile, -estimate_bin, -truth)
  
  print('Indexing samples according to order in data frame')
  IndexDF <- function(this_deeprvat_name){
    this_res = res %>% filter(model_name == this_deeprvat_name)
    this_res$index <- 1:nrow(this_res)
    return(this_res)
  }
  res = map_dfr(unique(res$model_name), IndexDF)
  stopifnot(length(unique(res %>% filter(index == 1) %>% pull(Y))) == 1 )
  this_plot_df = res %>% #filter((grepl(this_gene_list,model_name) | grepl("covariates",model_name))) %>%
    arrange(model_name, index)  %>%
    separate(col = model_name, into = c('burden', 'gene_list'), sep = '-', remove = FALSE) %>%
    rowwise() %>%
    mutate(model = ifelse(model_name == 'covariates-PRS', model_name, burden)) %>%
    ungroup() %>% 
    mutate(phenotype = this_phenotype) 
  
  unique(this_plot_df$model_name)
  
  # this_plot_df[['model']] = factor(this_plot_df[['model']], levels = scatter_plot_models_to_keep)
  
  plot_df_list[[this_phenotype]] = this_plot_df %>% mutate(phenotype = this_phenotype)
  
}

saveRDS(plot_df_list, file.path(out_dir, paste0('plot_df_list_', phenotype_suffix,'.Rds')))
