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

mkdir -p ./base
ln -s ../example_input_data_deeprvat/* ./base/


for rep in $(seq 0 $(($REPS - 1)) ); do
    mkdir -p ./sg_set_$rep
    ln -s ../example_input_data_deeprvat/* ./sg_set_$rep/
done

## ------------- CV SPLIT Version -------------------------------------------
# RUN CONFIG RULE - GENERATE CONFIG AND .PARQUET FILES ------------------------------------------
for rep in $(seq 0 $(($REPS - 1)) ); do
    now="$(date)"
    echo "Starting rep $rep: $now"
    cd ./sg_set_$rep
    snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile all_config -n
    cd ../
    now="$(date)"
    echo "Finished rep $rep: $now"
done

# RUN SEED GENE SET SELECTION ----------------------------------------------
# Run add_random_seed_genes.py to update seed_genes.parquet files for all Cv Splits in each rep dir with additional seed genes.
python $deeprvat_analysis_dir/add_random_seedgenes/add_random_seed_genes.py --reps $REPS --p_genes_random 0.2 --max_pval 0.5

# RUN CV ASSOCIATION TESTING PIPELINE ON EACH SEED GENE SET -------------------
for rep in $(seq 0 $(($REPS - 1)) ); do
    now="$(date)"
    echo "Starting rep $rep: $now"
    cd ./sg_set_$rep
    snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile --ri -n
    cd ../
    now="$(date)"
    echo "Finished rep $rep: $now"
done