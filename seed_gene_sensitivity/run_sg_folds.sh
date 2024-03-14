#!/usr/bin/env bash

# SETUP EXPERIMENT DIRECTORIES ----------------------------------------------
# Replace `[path_to_deeprvat]` with the path to your clone of the repository.
# Activate your virt-env before running script. i.e. $ mamba activate deeprvat

FOLDS=5
deeprvat_dir='/path/to/deeprvat/repo'
deeprvat_analysis_dir='/path/to/deeprvat_analysis/repo'

mkdir experiment_dir
cd ./experiment_dir

mkdir -p ./base
ln -s $deeprvat_dir/example/* ./base/

for fold in $(seq 0 $(($FOLDS - 1)) ); do
    mkdir -p ./fold_$fold
    ln -s $deeprvat_dir/example/* ./fold_$fold/
done

# RUN BASE EXP SEED GENE SELECTION ------------------------------------------
cd ./base
snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile all_config
cd ../

# RUN SEED GENE FOLD SELECTION ----------------------------------------------
# Run seed_gene_selection to yield seed_gene.parquet files in each FOLD dir
python $deeprvat_analysis_dir/seed_gene_sensitivity/seed_gene_selection.py --folds $FOLDS --downsample_percent 0.1 --min_keep_percent 0.5

# RUN ASSOCIATION TESTING PIPELINE ON EACH SEED GENE FOLD -------------------
for fold in $(seq 0 $(($FOLDS - 1)) ); do
    now="$(date)"
    echo "Starting Fold $fold: $now"
    cd ./fold_$fold
    snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile -n
    cd ../
    now="$(date)"
    echo "Finished Fold $fold: $now"
done