# Run phenotype prediction pipeline
The phenotype prediction pipeline uses cross-validation to attest generelizability of DeepRVAT-derived phenotype predictors. 
For this, the following steps are run on each CV fold. 
    1. Running Seed gene discory on training samples
    2. Training DeepRVAT on training samples using the seed genes from the fold
    3. Association testing using DeepRVAT gene impairments scores for the training samples 
    4. Fitting phenotype predictors using DeepRVAT gene impairment scores for the training samples for significantly associated genes from step 3
    5. Using the weights from the fitted predictors to predict the phenotype for test samples (using DeepRVAT gene impairment scores from the model trained in 3)
Afterwards, the predictions for test samples across all test folds are aggregated to evaluate the model. 


## Input data
An experiment directory with similar data input data as for the standard DeepRVAT pipeline has to be created (described [here](https://github.com/PMBio/deeprvat/)).
The annotations for seed gene discovery should be named `baseline_annotations.parquet`.
The annotations for DeepRVAT should be named `annotations.parquet`.
The `genes.parquet`, `phenotypes.parquet` (on the complete cohort), `variants.parquet` are the same as for the DeepRVAT pipeline.

### Genotype and phenotypes for each cross-validation fold
For running DeepRVAT using Cross-validation, the phenotype and genotype data has to be split by samples for each cv-fold. 
The split data has to be stored in `cv_data` in the experiment directory. 
The data in `cv_data` is the following: 
 `genotypes_train{x}.h5`, `genotypes_train{x}_phenotypes.parquet`, `genotypes_test{x}.h5`, `genotypes_test{x}_phenotypes.parquet`,
 where x = [0,1,2,3,4] for 5 CV folds.

`genotypes_train0.h5` for examples is a subset of `genotypes.parquet` from DeepRVAT, which only comprises the genotypes for the training samples from the first CV fold.


### PRS Scores
The common variant polygenic risk scores for all PGS ids listed in `../data/prs_pheno_map` have to be precomputed for the cohort. 
The resulting PRS have to be stored in `PRS.parquet` in the experiment directory. 
The column names of `PRS.parquet` have to be the PGS ids and the index (named `sample`) have to be the sample indices. 
`../data/prs_pheno_map` also has to be copied/linked to the experiment directory. 


## Run the pipeline
Copy `cv_deeprvat_training/config.yaml` and `phenotype_prediction_pipeline/config_eval.yaml` into the experiment directory.

`phenotype_prediction_pipeline/phenotype_prediction_pipeline.snakefile` can be used to run the enire pipeline. 
This will first run DeepRVAT and baseline methods in a cross-validation mode
using `cv_deeprvat_training/run_deeprvat_cv.snakefile`, which will use the config from `cv_deeprvat_training/config.yaml`

Then, the actual phenotype prediction models will be fitted using the following steps:

  1. Retrieve the burdens for all signficant genes from the results of the previous experimnet
  2. Prepare Regression model input data, comprising for all genes,
     the phenotype value, the PRS score, covariates (age, sex, genetic principal components) 
     and the gene scores (DeepRVAT or pLOF) for signficantly trait associated genes 
  3. Run the regression for logistic and linear models
  4. Prepare plotting data

When running the piplein, make sure the --use-conda flag is set


## Running the pipeline holding out 1 phenotype from DeepRVAT training

To generate the data for Figure 4 e&f,`phenotype_prediction_pipeline_hold_out_pheno.snakefile` has to be run. 
The pipeline assumes that the standard CV pipeline `run_deeprvat_cv.snakefile` has already been run. 
The `data_dir` in  `run_deeprvat_cv_hold_out_pheno.snakefile` and `phenotype_prediction_pipeline_hold_out_pheno.snakefile` has to be set to directory where `run_deeprvat_cv.snakefile` has been run. 
