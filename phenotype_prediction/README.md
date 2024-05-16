# Run phenotype prediction pipeline
This folder provides the code to run all analyses and generate all plots related to Figure 3 of the DeepRVAT manuscript. 

## Training/testing data sets
Training and evaluation of the phenotype prediction regression models is done on two disjoint data sets, restricting to unrelated Caucasian individuals. A total of 154,966 (from UKBB 200k WES) and 224,817 individuals (from UKBB 470k WES, not found in UKBB 200k WES) were used for training and evaluation, respectively.
The training data set is a subset of the unrelated Caucasian individuals from the UKBB 200k WES as used in the reference experiment, where samples with 3rd degree (or closer) relatives in the test set had been removed. 



## Input data
The experiment directory with the correct config has already been set-up `phenotype_prediction_exp`
The `protein_coding_genes.parquet`, `phenotypes.parquet` (on the complete cohort), `variants.parquet` are the same as for the DeepRVAT pipeline.

### Train/test sample files
A directory `train_test_samples` with the files `train_samples.pkl`, `test_samples.pkl`, comprising the sample ids for the training and testing fold (as a list of integers). 

### Genotype and phenotypes 
Genotypes (`genotypes.parquet`) and phenotypes (`phenotypes.parquet`) combined for all samples in  `train_test_samples` are required. 

### PRS Scores
The common variant polygenic risk scores for all PGS ids listed in `../data/prs_pheno_map` have to be precomputed for the cohort. 
The resulting PRS have to be stored in `PRS.parquet` in the experiment directory. 
The column names of `PRS.parquet` have to be the PGS ids and the index (named `sample`) have to be the sample ids. 
`../data/prs_pheno_map` also has to be copied/linked to the experiment directory. 

### DeepRVAT gene impairment scores
The folder `deeprvat_burdens` must have the `burdens.zarr` with the DeepRVAT gene impairment scores for all samples in `train_test_samples`. It also needs the metadata `sample_ids.zarr` and `genes.npy`, which provide the order of samples and genes in `burdens.zarr`. 
The gene impairment scores are computed with the pre-trained models from `../association_testing/paper_experiment`, which only used unrelated Caucasian samples from the  UKBB 200k WES cohort during training. 

### link to the reference experiment 
The lists of genes whose gene impairment scores will  be included in the rare variant phenotype predictors are retrieved from the reference experiment `../association_testing/paper_experiment`. 

### Alternative burden scores
The code to retrieve burden scores from single variant annotations is provided in `alternative_burden_computation`. It extracts the max (min for SIFT) annotation score for each sample and gene. 
`alternative_burden_computation/alternative_burdens.snakefile` has to be run before the `phenotype_prediction_pipeline.snakefile`. 
The `alternative_burden_computation` contain the the annotated variants, genotypes, and phenotypes for all samples
(`annotations.parquet`, `protein_coding_genes.parquet`, `genotypes.h5`, `phenotypes.parquet`, `variants.parquet`) (as [here](https://github.com/PMBio/deeprvat/example))


## Run the pipeline
From the experiment directory `phenotype_prediction_exp` run the `phenotype_prediction_pipeline.snakefile`, making sure that the `--use-conda` flag is set. 
This has to be run after `alternative_burden_computation/alternative_burdens.snakefile`!

The phenotype prediction pipeline will:
  1. Prepare Regression model input data, comprising for all genes,
     the phenotype value, the PRS score, covariates, 
     and the gene scores (DeepRVAT or or alternative burdens) for significantly trait associated genes. 
  2. Fit the regression models on the training data
  3. Evaluate the trained models on the testing data 
  4. Prepare plotting data


