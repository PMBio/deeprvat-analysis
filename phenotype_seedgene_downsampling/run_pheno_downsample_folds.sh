#!/usr/bin/env bash

# SETUP EXPERIMENT DIRECTORIES --------------------------------------------
# Prepare experiment directories with necessary input files for DeepRVAT running 
# Replace `[path_to_deeprvat]` with the path to your clone of the repository.
# Activate your virt-env before running script. i.e. $ mamba activate deeprvat

FOLDS=5
deeprvat_dir='/path/to/deeprvat/repo'
deeprvat_analysis_dir='/path/to/deeprvat_analysis/repo'

mkdir experiment_dir
cd ./experiment_dir

mkdir -p ./base
scp $deeprvat_dir/example/config.yaml ./base/config.yaml

for fold in $(seq 0 $(($FOLDS - 1)) ); do
    mkdir -p ./set_$fold
    ln -s $deeprvat_dir/example/* ./set_$fold/
done

# RUN PHENOTYPE SELECTION ----------------------------------------------
# Run phenotype_selection from experiment directory to yield new config.yaml files in each FOLD dir
python $deeprvat_analysis_dir/phenotype_sensitivity/phenotype_selection.py --folds $FOLDS --downsample_percent 0.1

# RUN CV ASSOCIATION TESTING PIPELINE ON EACH FOLD -----------------------
for fold in $(seq 0 $(($FOLDS - 1)) ); do
    now="$(date)"
    echo "Starting Fold $fold: $now"
    cd ./set_$fold
    snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile -n
    cd ../
    now="$(date)"
    echo "Finished Fold $fold: $now"
done