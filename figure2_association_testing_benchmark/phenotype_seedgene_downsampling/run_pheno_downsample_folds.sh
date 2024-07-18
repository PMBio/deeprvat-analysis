#!/usr/bin/env bash

# SETUP EXPERIMENT DIRECTORIES --------------------------------------------
# Prepare experiment directories with necessary input files for DeepRVAT running 
# Replace `[path_to_deeprvat]` with the path to your clone of the repository.
# Activate your virt-env before running script. i.e. $ mamba activate deeprvat

REPS=5
deeprvat_dir='/path/to/deeprvat/repo'
deeprvat_analysis_dir='/path/to/deeprvat_analysis/repo'

mkdir experiment_dir
cd ./experiment_dir

mkdir -p ./base
scp ../example_input_data_deeprvat/config.yaml ./base/config.yaml

for rep in $(seq 0 $(($REPS - 1)) ); do
    mkdir -p ./phenotype_ds/rep_$rep
    ln -s ../example_input_data_deeprvat/* ./phenotype_ds/rep_$rep/
done

# RUN PHENOTYPE SELECTION ----------------------------------------------
# Run phenotype_selection from experiment directory to yield new config.yaml files in each rep dir
python $deeprvat_analysis_dir/phenotype_sensitivity/phenotype_selection.py --reps $REPS --exp_dir ./phenotype_ds --downsample_percent 0.1

# RUN CV ASSOCIATION TESTING PIPELINE ON EACH rep -----------------------
for rep in $(seq 0 $(($REPS - 1)) ); do
    now="$(date)"
    echo "Starting rep $rep: $now"
    cd ./phenotype_ds/rep_$rep
    snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile -n
    cd ../
    now="$(date)"
    echo "Finished rep $rep: $now"
done