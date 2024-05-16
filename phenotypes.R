library(tibble)
#### Phenotype Names ####

monti_staar_name_dict = c('new_phenotypes' = '_new_phenotypes',
                          'all_phenotypes' =  '_all_phenotypes',
                          'training_phenotypes' = '_training_phenotypes')


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
                      "Red blood cell erythrocyte count" = "Erythrocyte count",
                      "Cholesterol statin corrected" = "Cholesterol",
                      "LDL direct statin corrected" = "LDL direct")

phenotype_renamer_rev = names(phenotype_renamer)
names(phenotype_renamer_rev) = phenotype_renamer

old_quant_phenotypes <- c(
  "Apolipoprotein_A",
  "Apolipoprotein_B",
  "Calcium",
  "Cholesterol_statin_corrected",
  "HDL_cholesterol",
  "IGF_1",
  "LDL_direct_statin_corrected",
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
  "Urate"
  )


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
    "Supraventricular tachycardia",
    "Tricuspid valve disease",
    "Ventricular arrhythmia",
    "Peripheral vascular disease",
    "Mitral valve disease",
    "Ischemic stroke",
    "Hypertrophic Cardiomyopathy",
    "Implantable cardioverter defibrillator",
    "Heart failure",
    "Dilated Cardiomyopathy",
    "Congenital heart disease",
    "Cardiac arrest",
    "Cardiac surgery",
    "Bradyarrhythmia",
    "Hypertension",
    "Hypercholesterolemia",
    "Coronary Artery Disease",
    "Atrial fibrillation",
    "Myocardial infarction",
    "Stroke",
    "Aortic valve disease",
    "AV or bundle branch block"
  ),
  "Respiratory" = c(
    "Sleep apnea",
    "Asthma",
    "Chronic obstructive pulmonary disease",
    "Pneumonia"
  ),
  "Gastrointestinal" = c(
    "Pancreatitis",
    "Inflammatory bowel disease",
    "Gastroesophageal reflux disease",
    "Diverticular disease",
    "Cholelithiasis",
    "Irritable bowel syndrome"
  ),
  "Neurological/Mental" = c( #Neurological and Mental Health
    "Parkinsons disease",
    "Epilepsy",
    "Anxiety",
    "Depression",
    "Migraine",
    "Bipolar Disorder"
  ),
  "Endocrine" = c(
    "Hyperthyroidism",
    "Diabetes Type 1",
    "Hypothyroidism",
    "Diabetes Type 2"
  ),
  "Dermatological" = c(
    "Skin cancer",
    "Dermatitis"
  ),
  "Ophthalmic" = c(
    "Glaucoma",
    "Cataract"
  ),
  "Allergies/Autoimmune" = c(
    "Psoriasis",
    "Multiple sclerosis",
    "Allergic rhinitis"
  ),
  "Musculoskeletal" = c(
    "Sciatica",
    "Rheumatoid arthritis",
    "Osteoarthritis",
    "Back pain",
    "Osteoporosis"
  ),
  "Cancer" = c(
    "Prostate cancer",
    "Lung cancer",
    "Colorectal cancer",
    "Cervical cancer",
    "Bladder cancer",
    "Breast cancer"
  ),
  "Venous" = c(
    "Venous thromboembolism"
  ),
  "Renal" = c(
    "Gout",
    "Chronic kidney disease"
  )
)

binary_trait_grouping_tibble <- enframe(binary_trait_grouping_list, name = "binary_trait_grouping", value = "phenotype") %>%
  unnest(cols = phenotype)

binary_phenotypes_imbalanced = c("Jurgens_Anxiety","Jurgens_Aortic_valve_disease","Jurgens_AV_or_bundle_branch_block","Jurgens_Bipolar_Disorder","Jurgens_Bladder_cancer","Jurgens_Bradyarrhythmia","Jurgens_Cardiac_arrest","Jurgens_Cardiac_surgery","Jurgens_Cervical_cancer","Jurgens_Chronic_kidney_disease","Jurgens_Colorectal_cancer","Jurgens_Congenital_heart_disease","Jurgens_Diabetes_Type_1","Jurgens_Dilated_Cardiomyopathy","Jurgens_Epilepsy","Jurgens_Glaucoma","Jurgens_Gout","Jurgens_Heart_failure","Jurgens_Hyperthyroidism","Jurgens_Hypertrophic_Cardiomyopathy","Jurgens_Implantable_cardioverter_defibrillator","Jurgens_Inflammatory_bowel_disease","Jurgens_Ischemic_stroke","Jurgens_Lung_cancer","Jurgens_Mitral_valve_disease","Jurgens_Multiple_sclerosis","Jurgens_Pancreatitis","Jurgens_Parkinsons_disease","Jurgens_Peripheral_vascular_disease","Jurgens_Prostate_cancer","Jurgens_Psoriasis","Jurgens_Rheumatoid_arthritis","Jurgens_Sciatica","Jurgens_Sleep_apnea","Jurgens_Supraventricular_tachycardia","Jurgens_Tricuspid_valve_disease","Jurgens_Ventricular_arrhythmia")

binary_phenotypes_balanced = c("Jurgens_Allergic_rhinitis","Jurgens_Asthma","Jurgens_Atrial_fibrillation","Jurgens_Back_pain","Jurgens_Breast_cancer","Jurgens_Cataract","Jurgens_Cholelithiasis","Jurgens_Chronic_obstructive_pulmonary_disease","Jurgens_Coronary_Artery_Disease","Jurgens_Depression","Jurgens_Dermatitis","Jurgens_Diabetes_Type_2","Jurgens_Diverticular_disease","Jurgens_Gastroesophageal_reflux_disease","Jurgens_Hypercholesterolemia","Jurgens_Hypertension","Jurgens_Hypothyroidism","Jurgens_Irritable_bowel_syndrome","Jurgens_Migraine","Jurgens_Myocardial_infarction","Jurgens_Osteoarthritis","Jurgens_Osteoporosis","Jurgens_Pneumonia","Jurgens_Skin_cancer","Jurgens_Stroke","Jurgens_Venous_thromboembolism")
binary_phenotypes = c(binary_phenotypes_imbalanced, binary_phenotypes_balanced)