
##################### Plotting.   ##################### 
colors <- c('#c6dbef', '#9ecae1', '#6baed6', '#3182bd', '#08519c', #'pink',
            "aquamarine4", '#E6AB02',
            '#7A2961')
methods <- c('Burden pLOF', 'Burden missense', 'SKAT pLOF', 'SKAT missense', 'Burden/SKAT combined', "Monti et al.", "STAAR", 'DeepRVAT')
main_methods = c('Burden/SKAT combined', "Monti et al.", "STAAR", 'DeepRVAT')
names(colors) = methods

binary_methods = c('SKAT missense',
               'SKAT pLOF', 
               'Burden missense', 
               'Burden pLOF',
               'Burden/SKAT combined',
               'DeepRVAT')


plot_dir = './plots'

font_size <- 10
title_font_size = 12
font_family <- "Helvetica"
plot_font <- element_text(size = font_size, family = font_family)

base_theme <- theme(
  text = element_text(size = font_size, family = font_family),
  axis.text = element_text(size = font_size, family = font_family),
  axis.title = element_text(size = font_size + 2, family = font_family),
  
  strip.text.x = element_text(size = font_size + 2, family = font_family),
  strip.text.y = element_text(size = font_size + 2, family = font_family),
  
  legend.title=element_blank(),
  plot.margin = margin(1, 1, 1, 1, "cm"),
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA, linewidth = 0), #transparent plot bg
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent', color=NA), #transparent legend panel
)
base_theme_wo_strip <- theme(
  text = element_text(size = font_size, family = font_family),
  axis.text = element_text(size = font_size, family = font_family),
  axis.title = element_text(size = font_size + 2, family = font_family),
  
  strip.text.x = element_blank(),
  strip.text.y = element_blank(),
  legend.title=element_blank(),
  plot.margin = margin(1, 1, 1, 1, "cm"),
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA, linewidth = 0), #transparent plot bg
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent', color=NA), #transparent legend panel
)


base_theme_wo_margin <- theme(
  text = element_text(size = font_size, family = font_family),
  axis.text = element_text(size = font_size, family = font_family),
  axis.title = element_text(size = font_size + 2, family = font_family),
  
  strip.text.x = element_text(size = font_size + 2, family = font_family),
  strip.text.y = element_text(size = font_size + 2, family = font_family),
  
  legend.title=element_blank(),
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA, linewidth = 0), #transparent plot bg
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent', color=NA), #transparent legend panel
)


checkResExistence = function(phenotypes, results_dir, results_dir_pattern = NA){
  all_files_exist = TRUE
  print(results_dir_pattern)
  for (p in phenotypes) {
    if(!is.na(results_dir_pattern)){
      this_res =  file.path(results_dir, p, results_dir_pattern, "all_results.parquet")
    }else{
      this_res = file.path(results_dir, p, "deeprvat", "eval", "all_results.parquet")
    }
    # Check if the file exists
    if (!file.exists(this_res)) {
      all_files_exist = FALSE
    }
  }
  return(all_files_exist)
}
############################################################ 

####################### Data functions ############################ 
loadDeepRVATResults = function(results_dir, phenotypes, phenotype_renamer, results_dir_pattern = NA){
  results_list <- list()
  print(paste('results_dir_pattern', results_dir_pattern))
  for (p in phenotypes) {
    print(p)
    if(!is.na(results_dir_pattern)){
      this_res =  read_parquet(file.path(results_dir, p, results_dir_pattern, "all_results.parquet"))
    }else{
      this_res =  read_parquet(file.path(results_dir, p, "deeprvat", "eval", "all_results.parquet"))
    }
    results_list[[p]] <- this_res 
  }
  results <- bind_rows(results_list) %>%
    rename("Trait" = "phenotype") %>% 
    mutate(Trait = ifelse(Trait %in% names(phenotype_renamer), phenotype_renamer[Trait], Trait)) %>%
    filter(Method %in% main_methods)

  
  significant <- results %>%
    filter(significant) %>%
    distinct(Trait, Method, gene, .keep_all = TRUE)

  # count sum of 'significant' with summarize instead of using '%>% count' since
  # count on significant genes only will loose trait/methods with zero counts
  counts_single = results %>% 
    select(Trait, Method, gene, significant) %>% 
    distinct() %>%
    group_by(Trait, Method) %>%
    summarise(n = sum(significant)) %>% 
    ungroup() %>% 
    mutate(is_single_trait = 'True')

  counts_all = counts_single %>%
    group_by(Method) %>%
    summarise(n = sum(n)) %>%
    ungroup() %>%
    mutate(Trait = "All traits", is_single_trait = "False") 

  counts = bind_rows(counts_single, counts_all)
  print('Trait/method combinations with zero counts')
  print(counts %>% filter(n == 0) )

  ##### Sanity check, can be removed in future #############################
  ###### the next part is just a sanity check to confirm that the counts as we count them now
  # (using group_by %>% summarise(n = sum(significant)) on results instead of count() on significant)
  # is the same as before for rows with non-zero n
  counts_old <- bind_rows(
    significant %>% mutate(is_single_trait = "True"),
    significant %>% mutate(Trait = "All traits", is_single_trait = "False")) %>%
    mutate(
      Method = factor(Method),
      Trait = factor(Trait)
    ) %>%
    count(Trait, Method, is_single_trait, .drop = FALSE)
  counts_old[is.na(counts_old$is_single_trait), "is_single_trait"] <- "True"
  df1 = counts %>% filter(n >0) 
  df2 = counts_old %>% filter(n >0) 
  stopifnot(all(dim(df1) == dim(df2)))
  merged_control = merge(df1, 
        df2, by = c("Trait", "Method", "is_single_trait"),
        suffixes = c("_df1", "_df2"))
  stopifnot(all(merged_control$n_df1 == merged_control$n_df2))
    
  # counts = counts %>% filter(Method != 'DeepRVAT wo baseline')
  ############# Sanity check end ####################################

  
  pheno_names_to_keep = as.character(unique(counts[['Trait']]))
  results = c('counts' = list(counts), 
              'results' = list(results), 
              'significant' = list(significant),
              'pheno_names_to_keep' = list(pheno_names_to_keep))
  return(results)
}

compute_lambda = function(grouped_df){
  lambda_gc <- grouped_df %>%
    arrange(pval) %>%
    summarise(
      observed_median_chi2 = median(pval, na.rm = TRUE),  # Step 1: Observed median chi-squared value
      expected_median_chi2 = qchisq(0.5, 1),             # Step 2: Expected median chi-squared value under the null hypothesis (df = 1)
      lambda_gc = observed_median_chi2 / expected_median_chi2  # Step 3: Genomic inflation factor (lambda_gc)
    )
  return(lambda_gc)
}

readReplicationData = function(replication_file_staar,
                               replication_file_monti,
                               replication_file_deeprvat){
  replication_staar = readRDS(replication_file_staar) %>% 
    mutate(Method = 'STAAR')
  replication_staar$Trait <- str_replace(replication_staar$Trait, "All Traits",  'All traits')
  replication_staar$Trait <- str_replace(replication_staar$Trait, "Mean platelet thrombocyte volume",  'MPTVS')
  
  replication_staar_all_traits =  replication_staar %>% 
    filter(Trait == 'All traits') %>% 
    rename(`Gene rank` = Rank, `Replicated genes` = Replicated)
  replication_staar_all_traits
  
  
  ## Monti
  replication_monti <- readRDS(replication_file_monti) %>%
    # filter(exp_name == "Monti-Rep Quantile-transformed") %>%
    mutate(Method = 'Monti et al.')
  
  replication_monti_all_traits <- replication_monti %>%
    filter(Trait == "All traits") %>%
    rename(`Gene rank` = Rank, `Replicated genes` = Replicated)
  replication_monti_all_traits
  
  ## DeepRVAT 
  replication <- read_parquet(replication_file_deeprvat)
  
  print('Differing trait names')
  print(symdiff(replication %>% distinct(Trait) %>% pull(Trait), replication_staar %>% distinct(Trait) %>% pull(Trait)))
  print(symdiff(replication %>% distinct(Trait) %>% pull(Trait), replication_monti %>% distinct(Trait) %>% pull(Trait)))
  
  replication <- replication %>% 
    filter(pheno_grouping == "all_phenotypes")
  replication %>% distinct(Trait)
  
  # Combining methods
  replication <- bind_rows(replication, replication_monti_all_traits, replication_staar_all_traits)
}



makeReplicationPlot = function(replication, max_rank, title = ''){

  replication$Significant <- recode(factor(replication$Significant, levels=c(TRUE, FALSE)), `TRUE` = "True", `FALSE` = "False")

  replication <- replication %>%
    mutate(Method = factor(Method, levels=methods)) %>%
    mutate(Trait = case_match(
      Trait,
      "MPTVS" ~ "MPTV",
      .default = Trait
    ))
  replication_plot <- ggplot(replication %>% filter(`Gene rank` < max_rank), aes(x=`Gene rank`, y=`Replicated genes`, color=Method)) +
    geom_abline(intercept=0, slope=1, color='lightgray', linewidth=0.5) +
    geom_step(linewidth=0.4, aes(linetype=Significant)) +
    scale_color_manual(values = colors) +
    labs(title = title) +
    # theme_cowplot() +
    theme_classic() +
    base_theme +
    theme(
      axis.text.x=element_text(angle=45, vjust=1, hjust=1),
      axis.text.y=element_text(angle=45),
      legend.position = "none"
    )
  
  max_significant <- replication %>% 
    filter(as.logical(Significant)) %>% 
    group_by(Method) %>% 
    summarize(`Gene rank` = max(`Gene rank`), `Replicated genes` = max(`Replicated genes`))
  max_significant
  replication_plot <- replication_plot + geom_point(data=max_significant, aes(x=`Gene rank`, y=`Replicated genes`, color=Method), size = 1)
  return(replication_plot)
}
