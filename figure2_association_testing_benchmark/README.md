# Instructions for reproducing results related to Figure. 2

### This directory provides code to reproduce all figures for the DeepRVAT association testing benchmark on unrelated caucasian indviduals from the 161k UKBB cohort (Figure 2 and related extended data figures).

All steps should be carried out using the deeprvat conda environment available in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).
Preprocess and annotate the UKBB WES data following the instructions in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).

## Run experiments

### Run Burden/SKAT baselines (seed gene discovery)
First run the seed gene discovery/burden & skat tests with their `config.yaml` files in `burden_skat` using the [seed gene discovery pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/seed_gene_discovery.snakefile). 

### Run DeepRVAT experiments
Once the Burden/SKAT baselines have is finsihed, DeepRVAT experiments can be run. This includes all directories starting with `deeprvat_`.
Each directory has to contain input data returned by the preprocessing and annotation pipelines (`genotypes.h5`, `variants.parquet`, `annotations.parquet`, `protein_coding_genes.parquet`, `phenotypes.parquet`). Also link these files into [example_input_data_deeprvat](example_input_data_deeprvat).

Then run the [DeepRVAT CV training/association testing pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/cv_training/cv_training_association_testing.snakefile) from all `deeprvat_` directories. 
One exception is the `deeprvat_synonymous` directory, where a symlink to the [pretrained_models](https://github.com/PMBio/deeprvat/tree/main/pretrained_models) from the DeepRVAT repository has to be created and association testing is subsequently run with the [pre-trained model association testing pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/association_testing_pretrained.snakefile).

To run the [phenotype and seed gene downsampling](phenotype_seedgene_downsampling) as well the experiments for assessing the  [robustness of the DeepRVAT to inclusion of null seed genes during training](add_random_seedgenes) run the shell scripts in the respective folders. 

For the conditional analysis, run the [conditional analysis pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/association_testing_control_for_common_variants.snakefile) from `deeprvat_main_exp` after the training/association testing has finished. 

### Run comparison methods
The experiments in `../comparison_methods/{monti/staar}/experiments` have to be run to get the results for Monti and STAAR. 

##  Plotting

### Prepare additional data required for plotting
To compute the replication results, run:
```
for EXP in deeprvat_main_exp deeprvat_linear_model deeprvat_plof_missense_anno
    do
    python compute_replication.py --out-dir $EXP $EXP
done
EXP=repeat_analysis && python compute_replication.py --analyze-all-repeats --out-dir $EXP $EXP
```

### Create the final plots

Run `figure2.Rmd`, `figure2_supp_downsampling_anno_removal.Rmd`,  `figure2_supp_random_seed_genes.Rmd`, `figure2_conditional_analysis.Rmd` to get all main and extended data figures for the association testing benchmark. 




