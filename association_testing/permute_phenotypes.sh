#!/usr/bin/env bash

python deeprvat_preprocess permute-phenotypes \
       --phenotype Apolipoprotein_A \
       --phenotype Apolipoprotein_B \
       --phenotype Calcium \
       --phenotype Cholesterol \
       --phenotype Red_blood_cell_erythrocyte_count \
       --phenotype HDL_cholesterol \
       --phenotype IGF_1 \
       --phenotype LDL_direct \
       --phenotype Lymphocyte_percentage \
       --phenotype Mean_platelet_thrombocyte_volume \
       --phenotype Mean_corpuscular_volume \
       --phenotype Mean_reticulocyte_volume \
       --phenotype Neutrophill_count \
       --phenotype Platelet_count \
       --phenotype Platelet_crit \
       --phenotype Platelet_distribution_width \
       --phenotype SHBG \
       --phenotype Standing_height \
       --phenotype Total_bilirubin \
       --phenotype Triglycerides \
       --phenotype Urate \
       paper_experiment/phenotypes.parquet \
       permutation_analysis/phenotypes.parquet
