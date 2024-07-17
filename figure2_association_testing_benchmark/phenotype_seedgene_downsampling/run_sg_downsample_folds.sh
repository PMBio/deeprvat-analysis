#!/usr/bin/env bash

# SETUP EXPERIMENT DIRECTORIES ----------------------------------------------
# Prepare experiment directories with necessary input files for DeepRVAT running 
# Replace `[path_to_deeprvat]` with the path to your clone of the repository.
# Activate your virt-env before running script. i.e. $ mamba activate deeprvat

REPS=5
deeprvat_dir='/path/to/deeprvat/repo'
deeprvat_analysis_dir='/path/to/deeprvat_analysis/repo'

mkdir experiment_dir
cd ./experiment_dir

mkdir -p ./base
ln -s $deeprvat_dir/example/* ./base/

for rep in $(seq 0 $(($REPS - 1)) ); do
    mkdir -p ./seedgene_ds/rep_$rep
    ln -s $deeprvat_dir/example/* ./seedgene_ds/rep_$rep/
done

# RUN BASE CV EXP SEED GENE SELECTION W/ CONFIG RULE ------------------------------------------
cd ./base
snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile all_config -n
cd ../

# RUN SEED GENE rep SELECTION ----------------------------------------------
# Run seed_gene_selection.py to yield seed_gene.parquet files in each rep dir
python $deeprvat_analysis_dir/phenotype_seedgene_downsampling/seed_gene_selection.py --reps $REPS --exp_dir ./seedgene_ds --downsample_percent 0.1 --min_keep_percent 0.5

# RUN CV ASSOCIATION TESTING PIPELINE ON EACH SEED GENE REP -------------------
for rep in $(seq 0 $(($REPS - 1)) ); do
    now="$(date)"
    echo "Starting rep $rep: $now"
    cd ./seedgene_ds/rep_$rep
    snakemake -j 1 --snakefile $deeprvat_dir/pipelines/cv_training/cv_training_association_testing.snakefile --ri -n
    cd ../
    now="$(date)"
    echo "Finished rep $rep: $now"
done