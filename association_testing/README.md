# Instructions for reproducing results related to Fig. 3 and associated supplementary figures

All steps should be carried out using the deeprvat conda environment available in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).

Preprocess and annotate the UKBB WES 200k release following the instructions in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).

## Run Burden/SKAT baselines (seed gene discovery)
First run the seed gene discovery/burden & skat tests with their `config.yaml` files in `burden_skat` and `burden_skat_binary` using the [seed gene discovery pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/seed_gene_discovery.snakefile). 

## Run DeepRVAT experiments for main figures

Once this is finsihed run the [DeepRVAT training/association testing pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/training_association_testing.snakefile), following the instructions in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/) in the `paper_experiment` folder (uses `paper_experiment/config.yaml` as the config file). 

To use the pre-trained DeepRVAT model from `paper_experiment` on quantitative and binary phenotypes that the model was not trained on, create a simlink `pretrained_models` pointing to `paper_experiment/models` in `paper_experiment` & `paper_experiment_binary`. 
Then run the [DeepRVAT association testing pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/association_testing_pretrained.snakefile) in `deeprvat_pretrained_quantitative` and  `deeprvat_pretrained_binary` 

Then run 
```
for x in $(\ls deeprvat_pretrained_quantitative | grep "^[A-Z]")
do
    ln -rs deeprvat_pretrained_quantitative/$x paper_experiment/$x
done
```
to link all the results for all quantitative phenotypes into one folder (required for the plotting markdown).

## Run DeepRVAT experiments for supplementary figures
 
For the subdirectories  `plof_missense_anno`, `linear_model`, `repeat_analysis`, with their association `config.yaml` files, run the [DeepRVAT training/association testing pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/training_association_testing.snakefile).

Before running the `permutation_analysis` experiment, run the script `permute_phenotypes.sh` to get a phenotypes dataframe with permuted phenotypes.

## Prepare additional data required for plotting
To compute the replication results, run:
```
for EXP in paper_experiment linear_model plof_missense_anno
    do
    python compute_replication.py --out-dir $EXP $EXP
done
EXP=repeat_analysis && python compute_replication.py --analyze-all-repeats --out-dir $EXP $EXP
```

## Get the paper figures

Use the notebooks `figure_3_main.Rmd` and `figure_3_supp.Rmd` to analyze the results. For this, the experiments in `../comparison_methods/{monti/staar}/experiments` also have to be run. 


