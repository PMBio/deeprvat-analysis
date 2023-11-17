#!/usr/bin/env bash

# SETUP EXPERIMENT DIRECTORIES ----------------------------------------------
# Replace `[path_to_deeprvat]` with the path to your clone of the repository.
# Activate your virt-env before running script. i.e. $ mamba activate deeprvat

FOLDS=5
path_to_deeprvat='/path/to/deeprvat/repo'

mkdir experiment_dir
cd ./experiment_dir

mkdir -p ./base
ln -s $path_to_deeprvat/example/* ./base/

for fold in $(seq 0 $(($FOLDS - 1)) ); do
    mkdir -p ./fold_$fold
    ln -s $path_to_deeprvat/example/* ./fold_$fold/
done

# RUN BASE EXP SEED GENE SELECTION ------------------------------------------
cd ./base
snakemake -j 1 --snakefile $path_to_deeprvat/pipelines/training_association_testing.snakefile all_config
cd ../

# RUN SEED GENE FOLD SELECTION ----------------------------------------------
# Run seed_gene_selection to yield seed_gene.parquet files in each FOLD dir
python ~/deeprvat-analysis-main/deeprvat-analysis/seed_gene_sensitivity/seed_gene_selection.py --folds $FOLDS --downsample_percent 0.1 --min_keep_percent 0.5

# RUN ASSOCIATION TESTING PIPELINE ON EACH SEED GENE FOLD -------------------
for fold in $(seq 0 $(($FOLDS - 1)) ); do
    cd ./fold_$fold
    snakemake -j 1 --snakefile $path_to_deeprvat/pipelines/training_association_testing.snakefile -n
    cd ../
done