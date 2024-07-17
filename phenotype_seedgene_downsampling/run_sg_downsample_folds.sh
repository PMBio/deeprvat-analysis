#!/usr/bin/env bash

# SETUP EXPERIMENT DIRECTORIES ----------------------------------------------
# Prepare experiment directories with necessary input files for DeepRVAT running 
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
    mkdir -p ./sg_set_$fold
    ln -s $deeprvat_dir/example/* ./sg_set_$fold/
done

# RUN BASE CV EXP SEED GENE SELECTION W/ CONFIG RULE ------------------------------------------
cd ./base
snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile all_config -n
cd ../

# RUN SEED GENE FOLD SELECTION ----------------------------------------------
# Run seed_gene_selection.py to yield seed_gene.parquet files in each FOLD dir
python $deeprvat_analysis_dir/phenotype_seedgene_downsampling/seed_gene_selection.py --folds $FOLDS --downsample_percent 0.1 --min_keep_percent 0.5

# RUN CV ASSOCIATION TESTING PIPELINE ON EACH SEED GENE FOLD -------------------
for fold in $(seq 0 $(($FOLDS - 1)) ); do
    now="$(date)"
    echo "Starting Fold $fold: $now"
    cd ./sg_set_$fold
    snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile --ri -n
    cd ../
    now="$(date)"
    echo "Finished Fold $fold: $now"
done