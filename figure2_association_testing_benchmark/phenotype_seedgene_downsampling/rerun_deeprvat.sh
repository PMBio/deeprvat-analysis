#!/usr/bin/env bash

# SETUP EXPERIMENT DIRECTORIES ----------------------------------------------
# Prepare experiment directories with necessary input files for DeepRVAT running
# Replace `[path_to_deeprvat]` with the path to your clone of the repository.
# Activate your virt-env before running script. i.e. $ mamba activate deeprvat

REPS=10
deeprvat_dir='/path/to/deeprvat/repo'
deeprvat_analysis_dir='/path/to/deeprvat_analysis/repo'

mkdir experiment_dir
cd ./experiment_dir



for rep in $(seq 0 $(($REPS - 1)) ); do
    mkdir -p ./rep_$rep
    ln -s ../example_input_data_deeprvat/* ./deeprvat_reruns/rep_$rep/
done


# RUN CV ASSOCIATION TESTING PIPELINE ON EACH REPLICATE -------------------
# This runs multiple runs of the "standard" DeepRVAT cv training + association testing (as in deeprvat-analysis/association_testing/deeprvat_main_exp) 
# to asses the random variability of DeepRVAT 

for rep in $(seq 0 $(($REPS - 1)) ); do
    now="$(date)"
    echo "Starting rep $rep: $now"
    cd ./deeprvat_reruns/rep_$rep
    snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile  -n
    cd ../
    now="$(date)"
    echo "Finished rep $rep: $now"
done