
#### Phenotype Names ####
binary_phenotypes = c(
  "Jurgens_Hypertension", 
  "Jurgens_Hypercholesterolemia", 
  "Jurgens_Osteoarthritis", 
  "Jurgens_Asthma",
  "Jurgens_Gastroesophageal_reflux_disease",
  "Jurgens_Diverticular_disease",
  "Jurgens_Depression", 
  "Jurgens_Cataract", 
  "Jurgens_Allergic_rhinitis", 
  "Jurgens_Hypothyroidism", 
  "Jurgens_Diabetes_Type_2", 
  "Jurgens_Skin_cancer",
  "Jurgens_Back_pain",
  "Jurgens_Cholelithiasis", 
  "Jurgens_Pneumonia",
  "Jurgens_Atrial_fibrillation", 
  "Jurgens_Coronary_Artery_Disease",
  "Jurgens_Migraine",
  "Jurgens_Dermatitis", 
  "Jurgens_Irritable_bowel_syndrome", 
  "Jurgens_Chronic_obstructive_pulmonary_disease", 
  "Jurgens_Venous_thromboembolism", 
  "Jurgens_Breast_cancer", 
  "Jurgens_Osteoporosis", 
  "Jurgens_Myocardial_infarction", 
  "Jurgens_Stroke" 
)

monti_staar_name_dict = c('new_phenotypes' = '_new_phenotypes',
                          'all_phenotypes' =  '_all_phenotypes',
                          'all' = '')


phenotype_renamer = c('WHR Body mass index BMI corrected' = 'WHR', 
                      'Forced expiratory volume in 1 second FEV1' = 'FEV1', 
                      'Glycated haemoglobin HbA1c' = 'HbA1c',
                      'Body mass index BMI' = 'BMI',
                      'Alkaline phosphatase' = 'ALP', 
                      'Gamma glutamyltransferase' = 'GGT',
                      'Whole body fat free mass' = 'Fat free mass',
                      "IGF 1" = "IGF-1",
                      "Mean platelet thrombocyte volume" = "MPTV",
                      "MPTVS" = "MPTV",
                      "Red blood cell erythrocyte count" = "Erythrocyte count")

phenotype_renamer_rev = names(phenotype_renamer)
names(phenotype_renamer_rev) = phenotype_renamer

old_quant_phenotypes <- c(
  "Apolipoprotein_A",
  "Apolipoprotein_B",
  "Calcium",
  "Cholesterol",
  "HDL_cholesterol",
  "IGF_1",
  "LDL_direct",
  "Lymphocyte_percentage",
  "Mean_corpuscular_volume",
  "Mean_platelet_thrombocyte_volume",
  "Mean_reticulocyte_volume",
  "Neutrophill_count",
  "Platelet_count",
  "Platelet_crit",
  "Platelet_distribution_width",
  "Red_blood_cell_erythrocyte_count",
  "SHBG",
  "Standing_height",
  "Total_bilirubin",
  "Triglycerides",
  "Urate")


new_quant_phenotypes = c(
  "Body_mass_index_BMI",
  "Glucose",
  "Vitamin_D",
  "Albumin",
  "Total_protein",
  "Cystatin_C",
  "Gamma_glutamyltransferase",
  "Alkaline_phosphatase",
  "Creatinine",
  "Whole_body_fat_free_mass", 
  "Forced_expiratory_volume_in_1_second_FEV1",
  "Glycated_haemoglobin_HbA1c",
  "WHR_Body_mass_index_BMI_corrected"
)

old_phenotype_names = gsub('_', ' ', old_quant_phenotypes)
old_phenotype_names = ifelse(old_phenotype_names %in% names(phenotype_renamer), phenotype_renamer[old_phenotype_names], old_phenotype_names)
new_phenotype_names = gsub('_', ' ', new_quant_phenotypes)
new_phenotype_names = ifelse(new_phenotype_names %in% names(phenotype_renamer), phenotype_renamer[new_phenotype_names], new_phenotype_names)

names(new_quant_phenotypes) = new_phenotype_names
names(old_quant_phenotypes) = old_phenotype_names


###
phenotype_old_new_dict <- c()
for (i in old_quant_phenotypes) {
  phenotype_old_new_dict[i] <- "Training Trait"
}

for (i in new_quant_phenotypes) {
  phenotype_old_new_dict[i] <- "New Trait"
}


new_names = names(phenotype_old_new_dict)
new_names = gsub('_', ' ', names(phenotype_old_new_dict))
matching_indices <- match(new_names, names(phenotype_renamer))
new_names[!is.na(matching_indices)] <- phenotype_renamer[matching_indices[!is.na(matching_indices)]]
names(phenotype_old_new_dict) = new_names

quant_phenotypes = c(old_quant_phenotypes, new_quant_phenotypes)


## Trait groupings 

quant_trait_grouping = c(
  "Bone and joint" = c("Vitamin D",  "ALP", "Calcium"),
  "Lipids" = c( "Apolipoprotein A","Apolipoprotein B",
                "Cholesterol","HDL cholesterol","Triglycerides", "LDL direct"),
  "Renal"  = c("Total protein", "Creatinine", "Cystatin C", "Urate"),
  "Liver" = c("Albumin" , "GGT", "Total bilirubin"),
  "Diabetes" = c("HbA1c", "Glucose"),
  "Physical measures" = c("BMI", "WHR", "Fat free mass", 'FEV1', "Standing height"),
  "Hormonal" = c("SHBG", "IGF-1"),
  "Blood count" = c( "Lymphocyte percentage","Mean corpuscular volume","MPTV", "Mean reticulocyte volume","Neutrophill count","Platelet count","Platelet crit","Platelet distribution width","Erythrocyte count")
)

binary_trait_grouping_list <- list(
  "Cardiovascular" = c(
    "Hypertension",
    "Hypercholesterolemia",
    "Coronary Artery Disease",
    "Atrial fibrillation",
    "Myocardial infarction",
    "Stroke"
  ),
  "Respiratory" = c(
    "Asthma",
    "Chronic obstructive pulmonary disease",
    "Pneumonia"
  ),
  "Gastrointestinal" = c(
    "Gastroesophageal reflux disease",
    "Diverticular disease",
    "Cholelithiasis",
    "Irritable bowel syndrome"
  ),
  "Neurological and Mental Health" = c(
    "Depression",
    "Migraine"
  ),
  "Endocrine" = c(
    "Hypothyroidism",
    "Diabetes Type 2"
  ),
  "Dermatological" = c(
    "Skin cancer",
    "Dermatitis"
  ),
  "Ophthalmic" = c(
    "Cataract"
  ),
  "Allergies and Immune" = c(
    "Allergic rhinitis"
  ),
  "Musculoskeletal" = c(
    "Osteoarthritis",
    "Back pain",
    "Osteoporosis"
  ),
  "Cancer" = c(
    "Breast cancer"
  ),
  "Venous" = c(
    "Venous thromboembolism"
  )
)
binary_trait_grouping_tibble <- enframe(binary_trait_grouping_list, name = "binary_trait_grouping", value = "Trait") %>%
  unnest(cols = Trait)


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
############################################################ 


####################### Data functions ############################ 
loadDeepRVATResults = function(results_dir, phenotypes, phenotype_renamer){
  results_list <- list()
  for (p in phenotypes) {
    print(p)
    this_res =  read_parquet(file.path(results_dir, p, "deeprvat", "eval", "all_results.parquet"))
    results_list[[p]] <- this_res 
  }
  results <- bind_rows(results_list) %>%
    rename("Trait" = "phenotype") %>% 
    mutate(Trait = ifelse(Trait %in% names(phenotype_renamer), phenotype_renamer[Trait], Trait)) %>%
    filter(Method %in% main_methods)

  
  significant <- results %>%
    filter(significant) %>%
    distinct(Trait, Method, gene, .keep_all = TRUE)
  
  counts <- bind_rows(
    significant %>% mutate(is_single_trait = "True"),
    significant %>% mutate(Trait = "All traits", is_single_trait = "False")) %>%
    mutate(
      Method = factor(Method),
      Trait = factor(Trait)
    ) %>%
    count(Trait, Method, is_single_trait, .drop = FALSE)
  counts[is.na(counts$is_single_trait), "is_single_trait"] <- "True"
  # counts = counts %>% filter(Method != 'DeepRVAT wo baseline')
  
  pheno_names_to_keep = as.character(unique(counts[['Trait']]))
  results = c('counts' = list(counts), 
              'results' = list(results), 
              'significant' = list(significant),
              'pheno_names_to_keep' = list(pheno_names_to_keep))
  return(results)
}

compute_lambda = function(grouped_df){
  lambda_gc <- grouped_df %>%
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
  replication_plot <- replication_plot + geom_point(data=max_significant, aes(x=`Gene rank`, y=`Replicated genes`, color=Method))
  return(replication_plot)
}
