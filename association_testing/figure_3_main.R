library(ggplot2)
library(arrow)
library(dplyr)
library(stringr)
library(tidyr)
library(ggpubr)
library(cowplot)
library(tibble)

# code_dir = '/home/e400p/deeprvat_public/deeprvat-analysis/association_testing'
code_dir = snakemake@params[["code_dir"]]
source(file.path(code_dir, 'rvat_figure_utils.R'))
source(file.path(code_dir, 'rvat_figure_paths.R'))



## Define experiment directories
# replication_file_deeprvat = '~/ukbb/experiments/experiments_eva/replication_dataset_bugfix/deeprvat_alpha_mvc0/replication.parquet'
# results_dir = '/home/e400p/ukbb/experiments/revision_2/deeprvat_cv/h0_seed_genes_10_5_fold_standard_config'
# results_dir = "~/ukbb/experiments/dataset_bugfix/deeprvat_alpha_mvc0/"
replication_file_deeprvat = snakemake@input[["replication"]]
results_dir = snakemake@params[["results_dir"]]
out_file = snakemake@output[[1]][1]
print(out_file)
results_dir_pattern = snakemake@params[["results_dir_pattern"]]
results_dir_pattern = ifelse(is.null(results_dir_pattern), NA, results_dir_pattern)

print(sprintf('results dir pattern %s', results_dir_pattern))

print(sprintf('DeepRVAT quant result dir: %s: replication file deeprvat %s', results_dir, replication_file_deeprvat))

cat('Other result dirs: ', monti_res_binary_dir,
    monti_quant_dir,
    staar_quant_dir,
    staar_res_binary_dir, 
    results_deeprvat_binary_dir,
    timing_data_dir, sep = '\n')


if(checkResExistence(quant_phenotypes, results_dir, results_dir_pattern)){
  pheno_choice = 'all_phenotypes'
}else if (checkResExistence(old_quant_phenotypes, results_dir, results_dir_pattern)){
  pheno_choice = 'training_phenotypes'
}else if (checkResExistence(new_quant_phenotypes, results_dir, results_dir_pattern)){
  pheno_choice = 'new_phenotypes'
} else{
  stop('files for all options of pheno_choice are missing')
}
print(sprintf('Analyzing for pheno choice %s', pheno_choice))


quant_pheno_dict = c('all_phenotypes' = list(quant_phenotypes), 
                     'training_phenotypes' = list(old_quant_phenotypes),
                     'new_phenotypes' = list(new_quant_phenotypes))

quant_phenotypes_to_analyse = quant_pheno_dict[[pheno_choice]]

## DeepRVAT


all_deeprvat_res = loadDeepRVATResults(results_dir = results_dir, 
  phenotypes = quant_phenotypes_to_analyse, 
  phenotype_renamer = phenotype_renamer,
  results_dir_pattern = results_dir_pattern)
counts_deeprvat_quant = all_deeprvat_res[['counts']]
significant_deeprvat_quant = all_deeprvat_res[['significant']]
results_deeprvat_quant = all_deeprvat_res[['results']]
pheno_names_to_keep = all_deeprvat_res[['pheno_names_to_keep']]



## Load Monti et al. results


counts_monti_quant = readRDS(file.path(monti_quant_dir, paste0("monti_counts", monti_staar_name_dict[pheno_choice], ".Rds"))) %>%
  filter(Trait != 'WHR') %>%
  mutate(is_single_trait = ifelse(Trait == 'All traits', 'False', 'True')) %>%
  rename(n = discoveries) %>%
  mutate(Method = 'Monti et al.') %>%
  mutate(Trait = ifelse(Trait %in% names(phenotype_renamer), phenotype_renamer[Trait], Trait)) %>%
  filter(Trait %in% pheno_names_to_keep)

monti_res_quant = readRDS(file.path(monti_quant_dir, paste0("monti_results", monti_staar_name_dict[pheno_choice], ".Rds")))  %>% 
  filter(Trait != 'WHR') %>%
  mutate(Method = 'Monti et al.') %>%
  mutate(Trait = ifelse(Trait %in% names(phenotype_renamer), phenotype_renamer[Trait], Trait)) %>% 
  filter(Trait %in% pheno_names_to_keep)
monti_res_quant %>% distinct(Trait)


assertthat::assert_that(setequal(counts_deeprvat_quant %>% distinct(Trait) %>% pull(Trait), counts_monti_quant %>% distinct(Trait) %>% pull(Trait)))
assertthat::assert_that(setequal(results_deeprvat_quant %>% distinct(Trait) %>% pull(Trait), monti_res_quant %>% distinct(Trait) %>% pull(Trait)))

all_traits_count = sum(counts_monti_quant %>% filter(is_single_trait == 'True') %>% pull(n))
counts_monti_quant = counts_monti_quant %>% mutate(n = ifelse(Trait == 'All traits',all_traits_count, n))


## Load STAAR results



counts_staar_quant = readRDS(file.path(staar_quant_dir, paste0("sig_counts_staar", monti_staar_name_dict[pheno_choice], ".Rds"))) %>% 
  filter(Trait != 'WHR') %>%
  mutate(Method = 'STAAR')  %>%
  rename(n = discoveries) %>%
  mutate(is_single_trait = ifelse(Trait == 'All traits', 'False', 'True')) %>%
  mutate(Trait = ifelse(Trait %in% names(phenotype_renamer), phenotype_renamer[Trait], Trait)) %>%
  filter(Trait %in% pheno_names_to_keep)

all_traits_count = sum(counts_staar_quant %>% filter(is_single_trait == 'True') %>% pull(n))
counts_staar_quant = counts_staar_quant %>% mutate(n = ifelse(Trait == 'All traits',all_traits_count, n))


staar_res_quant = readRDS(file.path(staar_quant_dir, paste0("staar_corrected", monti_staar_name_dict[pheno_choice], ".Rds")))  %>% 
  filter(Trait != 'WHR') %>%
  mutate(Method = 'STAAR') %>%
  mutate(Trait = ifelse(Trait %in% names(phenotype_renamer), phenotype_renamer[Trait], Trait)) %>% 
  filter(Trait %in% pheno_names_to_keep)  %>%
  drop_na(pval) %>%
  mutate(pval = as.numeric(pval))



assertthat::assert_that(setequal(counts_deeprvat_quant %>% distinct(Trait) %>% pull(Trait), counts_staar_quant %>% distinct(Trait) %>% pull(Trait)))
assertthat::assert_that(setequal(results_deeprvat_quant %>% distinct(Trait) %>% pull(Trait), staar_res_quant %>% distinct(Trait) %>% pull(Trait)))




## Combine quantiative counts results

 
all_counts_quant =  rbind(counts_deeprvat_quant, counts_monti_quant, counts_staar_quant) %>%
  mutate(Method = factor(Method, levels=methods)) %>%
  mutate(Trait = as.character(Trait)) %>%
  mutate(Trait = recode(Trait, !!!phenotype_renamer))  %>%  
  mutate(pheno_group = case_when(
    Trait %in% names(old_quant_phenotypes) ~ 'old',
    Trait %in% names(new_quant_phenotypes) ~ 'new',
    TRUE ~ 'all'
  )) 

all_counts_quant %>% distinct(Trait, pheno_group)
all_counts_quant %>% filter(Trait == "All traits")
all_counts_quant %>% filter(Trait != "All traits")
levels(all_counts_quant$Method)
all_counts_quant



all_counts_quant %>% distinct(pheno_group)
all_counts_quant %>% filter(pheno_group == 'Training Trait') %>% 
  filter(Method != 'STAAR')


## group quantiative counts


flattened_quant_trait_grouping <- stack(quant_trait_grouping) %>%
  mutate(ind = gsub('\\d+', '', ind))
quant_trait_grouping = flattened_quant_trait_grouping [['ind']]
names(quant_trait_grouping) = flattened_quant_trait_grouping[['values']]
quant_trait_grouping

setdiff(all_counts_quant %>% distinct(Trait) %>% pull(Trait), names(quant_trait_grouping))
all_counts_quant = all_counts_quant %>% mutate(trait_group = quant_trait_grouping[as.character(Trait)]) %>%
  mutate(pheno_group = phenotype_old_new_dict[as.character(Trait)]) %>%
  mutate(pheno_group =  factor(pheno_group, levels = c('Training Trait', 'New Trait')))

grouped_counts_by_pheno_group = all_counts_quant %>% filter(is_single_trait == 'True') %>%
  group_by(Method, trait_group, pheno_group) %>%
  summarise(n = sum(n)) %>%
  ungroup()
grouped_counts_all_pheno_group = all_counts_quant %>% filter(is_single_trait == 'True') %>%
  group_by(Method, pheno_group) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  mutate(trait_group = 'All Traits')
grouped_counts = rbind(grouped_counts_by_pheno_group, grouped_counts_all_pheno_group)  
grouped_counts %>% filter(pheno_group == 'All Traits')

trait_counts = all_counts_quant %>% filter(is_single_trait == 'True') %>%
  select(pheno_group, Trait, trait_group) %>%
  distinct() %>%
  group_by(trait_group, pheno_group) %>%
  summarise(n_phenotypes = n()) %>%
  mutate(trait_with_count = paste0(trait_group, ' (',n_phenotypes, ')' ))

trait_counts = rbind(trait_counts, all_counts_quant %>% filter(is_single_trait == 'True') %>%
                       select(pheno_group, Trait) %>%
                       distinct() %>%
                       group_by(pheno_group) %>%
                       summarise(n_phenotypes = n()) %>% 
                       mutate(trait_group = 'All Traits')) %>%
  mutate(trait_with_count = paste0(trait_group, ' (',n_phenotypes, ')' ))

trait_counts

grouped_counts = grouped_counts %>% left_join(trait_counts)
grouped_counts = grouped_counts %>% mutate(is_single_trait = ifelse(trait_group  == 'All Traits', FALSE, TRUE)) 


# Figure 3 plots


plot_list = list()


trait_title_dict = c('Training Trait' = 'Traits used in gene impairment module training', 'New Trait' = 'Traits not used in gene impairment module training')

t_groups_dict = c('all_phenotypes' =  list(c('Training Trait', 'New Trait')), 
                     'training_phenotypes' = c('Training Trait'),
                     'new_phenotypes' = c('New Trait'))
print('plotting')
for (t_group in t_groups_dict[[pheno_choice]]){
  p = ggplot(
    grouped_counts %>% filter(pheno_group == t_group & Method %in% main_methods),
    aes(
      x=trait_with_count,
      y=n,
      fill=Method,
    )
  ) + 
    geom_col(position='dodge', width=0.9, color='darkgray', linewidth=0.1) +
    labs(y = "Significant gene-trait associations", title = trait_title_dict[[t_group]]) +
    theme_classic() +
    # facet_grid(.~is_single_trait, scales='free', space = 'free') + #,  space = "free_x") +
    ggforce::facet_row(vars(is_single_trait), scales = 'free', space = 'free') +
    # theme_cowplot() +
    base_theme_wo_strip +
    theme(
      axis.text.x=element_text(angle=45, vjust=1, hjust=1),
      plot.title = element_text(hjust = 0.5, size = title_font_size),
      axis.title.x=element_blank(),
      legend.key.size=unit(0.5, "line"),
      legend.text = element_text(size = 12),
      panel.spacing=unit(0.75, "cm"),
      strip.background = element_blank(),
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) 
  print(p)
  plot_list[[paste('quant_discoveries', gsub(' ', '_', t_group), sep = "_")]] = p
}
grouped_counts %>% filter(trait_group == 'All Traits')



# Replication training traits
# Here we only check the replication on the training phenotypes

pheno_suffix = '_training_phenotypes'

replication_file_staar = file.path(staar_quant_dir, paste0('replication_staar', pheno_suffix, '.Rds'))
replication_file_monti = file.path(monti_quant_dir, paste0('replication_monti', pheno_suffix, '.Rds'))
replication_data = readReplicationData(replication_file_staar, replication_file_monti, replication_file_deeprvat)
replication_plot = makeReplicationPlot(replication_data %>% filter(Method %in% main_methods), max_rank = 1000)
replication_plot

plot_list[['replication_training']] = replication_plot

print('plotting')
if (pheno_choice == 'all_phenotypes'){
  ggsave(out_file,
         ggarrange(plot_list[['quant_discoveries_Training_Trait']], plot_list[['replication_training']], plot_list[['quant_discoveries_New_Trait']],
                   nrow = 2, ncol =2 , common.legend = TRUE),
         width = 973/90, height = 949/90)
}else{
  ggsave(out_file,
         ggarrange(plot_list[[paste('quant_discoveries', gsub(' ', '_', t_groups_dict[[pheno_choice]]), sep = '_')]], plot_list[['replication_training']],
                   nrow = 2, ncol =2 , common.legend = TRUE),
         width = 973/90, height = 949/90)
}

print('done')




