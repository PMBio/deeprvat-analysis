### Instructions for reproducing results related to Fig. 3 and associated supplementary figures

All steps should be carried out using the deeprvat conda environment available in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).

Preprocess and annotate the UKBB WES 200k release following the instructions in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).


First run the seed gene discovery/burden & skat tests with their `config.yaml` files in `burden_skat` and `burden_skat_binary` using the [seed gene discovery pipeline] https://github.com/PMBio/deeprvat/blob/main/pipelines/seed_gene_discovery.snakefile. 

Once this is finsihed, for the subdirectories `paper_experiment`, `plof_missense_anno`, `linear_model`, `repeat_analysis`, with their association `config.yaml` files, run the [DeepRVAT training/association testing pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/training_association_testing.snakefile), following the instructions in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/). 

Before running the `permutation_analysis` experiment, run the script `permute_phenotypes.sh` to get a phenotypes dataframe with permuted phenotypes.

To use the pre-trained DeepRVAT model from `paper_experiment` on quantitative and binary phenotypes that the model was not trained on, create a simlink `pretrained_models` pointing to `paper_experiment/models` in `paper_experiment` & `paper_experiment_binary`. 
Then run this [pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/association_testing_pretrained.snakefile) in `paper_experiment` using `config_pretrained_quant_phenotypes.yaml` and in `paper_experiment_binary`.


Compute the replication results:
```
for EXP in paper_experiment linear_model plof_missense_anno
    do
    python compute_replication.py --out-dir $EXP $EXP
done
EXP=repeat_analysis && python compute_replication.py --analyze-all-repeats --out-dir $EXP $EXP
```

Use the notebooks `figure_3_main.Rmd` and `figure_3_supp.Rmd` to analyze the results. For this, the experiments in `../comparison_methods` also have to be run. 

