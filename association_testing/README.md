### Instructions for reproducing results related to Fig. 3 and associated supplementary figures

All steps should be carried out using the deeprvat conda environment available in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).

Preprocess and annotate the UKBB WES 200k releasefollowing the instructions in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/).

For the three subdirectories with their association `config.yaml` files, run the seed gene discovery and DeepRVAT training/association testing pipelines, following the instructions in the [main DeepRVAT repository](https://github.com/PMBio/deeprvat/). Before running the `permutation_analysis` experiment, run the script `permute_phenotypes.sh` to get a phenotypes dataframe with permuted phenotypes.

Compute the replication results:
```
for EXP in paper_experiment linear_model plof_missense_anno permutation_analysis
    do
    python compute_replication.py --out-dir $EXP $EXP
done
```

Use the notebook `rvat_figures.Rmd` to analyze the results.
