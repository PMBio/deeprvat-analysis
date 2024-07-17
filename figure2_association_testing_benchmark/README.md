# Instructions for reproducing results related to association testing figures from the DeepRVAT paper (Fig. 2 and Fig.4) and associated supplementary figures

All steps should be carried out using the deeprvat conda environment available in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).

Preprocess and annotate the UKBB WES data following the instructions in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).

## Run Burden/SKAT baselines (seed gene discovery)
First run the seed gene discovery/burden & skat tests with their `config.yaml` files in `burden_skat` using the [seed gene discovery pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/seed_gene_discovery.snakefile). 

## Run DeepRVAT experiments for main figures

Once this is finsihed run the [DeepRVAT CV training/association testing pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/cv_training/cv_training_association_testing.snakefile), following the instructions in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/) in the `deeprvat_main_exp` folder (uses `deeprvat_main_exp/config.yaml` as the config file). 

Consider the instructions for the [CV training](https://deeprvat.readthedocs.io/en/latest/deeprvat.html#training-and-association-testing-using-cross-validation) to run the experiment. 


## Run DeepRVAT experiments for supplementary figures
 
For the subdirectories  `deeprvat_plof_missense_anno`, `deeprvat_linear_model`, with their association `config.yaml` files, run the  [DeepRVAT CV training/association testing pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/cv_training/cv_training_association_testing.snakefile).


## Prepare additional data required for plotting
To compute the replication results, run:
```
for EXP in deeprvat_main_exp linear_model plof_missense_anno
    do
    python compute_replication.py --out-dir $EXP $EXP
done
EXP=repeat_analysis && python compute_replication.py --analyze-all-repeats --out-dir $EXP $EXP
```
## Experiments for Figure 4
The DeepRVAT model used to compute DeepRVAT gene impairment scores for all analyses in Figure 4 is the one trained in `deeprvat_main_exp`. 
Follow the instructions from the main DeepRVAT instructions for using [DeepRVAT with REGENIE](https://deeprvat.readthedocs.io/en/latest/deeprvat.html#running-the-association-testing-pipeline-with-regenie) with precomputed gene impairment scores for quantiative and binary traits and on diverse subsets of samples. 
For running the default REGENIE RVAT test (burden/SKAT) see the [REGENIE documentation](https://rgcgithub.github.io/regenie/). 

## Get the paper figures

Use the makrdowns in this directory to generate all figures. For this, the experiments in `../comparison_methods/{monti/staar}/experiments` also have to be run. 
