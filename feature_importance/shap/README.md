## SHAP experiments

### Generate supplementary figures

You can use [explainer.ipynb](https://github.com/PMBio/deeprvat-analysis/blob/main/feature_importance/mutagenesis/mutagenesis.ipynb) to generate Supp  Figure 3.5 (b) and (c) using [shap_relative_per_sampling.pkl](https://github.com/PMBio/deeprvat-analysis/blob/main/feature_importance/shap/shap_relative_per_sampling.pkl)  and [shap_score_per_sampling.pkl](https://github.com/PMBio/deeprvat-analysis/blob/main/feature_importance/shap/shap_score_per_sampling.pkl) files. 

### Running snakemake pipeline

Use `explain.snakefile` to run SHAP experiments across different samplings to generate raw files which are used to obtain the aggreagted pickle files above. 

To be able to run this pipeline, you need to replace the following variables in `explain.snakefile` file:
* `[path_to_deeprvat]` with the path to your clone of [deeprvat](https://github.com/PMBio/deeprvat/) repository  
* `[path_to_deeprvat-analysis]`  with the path to your clone of this repository  
* `[path_to_inputs]` with the path that contain the outputs of complete model training where for each phenotype following four files are generated:
  * `input_tensor_file = [path_to_inputs]/{pheno}/deeprvat/input_tensor.zarr`
  * `covariates_file = [path_to_inputs]/{pheno}/deeprvat/covariates.zarr`
  * `y_file = [path_to_inputs]/{pheno}/deeprvat/y.zarr`
  * `training_gene_file = [path_to_inputs]/{pheno}/deeprvat/seed_genes.parquet`




### Running [SHAP DeepExplainer](https://github.com/slundberg/shap)

Make sure SHAP package is installed in your [deeprvat_environment].

```bash
	conda activate [deeprvat_environment]
	conda install -c conda-forge shap==0.42.1 
```

SHAP-explainer pipeline uses a modified version of [models.py](https://github.com/PMBio/deeprvat/blob/master/deeprvat/deeprvat/models.py) file to accommadate the design of SHAP Deep Explainer. Therefore in order to run SHAP you need to replace the original `model.py` with the one provided in this directory. 


```bash
	conda activate [deeprvat_environment]
	cd [path_to_deeprvat_analysis]/feature_importance/shap/
	mv [path_to_deeprvat]/deeprvat/deeprvat/deeprvat/models.py [path_to_deeprvat]/deeprvat/deeprvat/deeprvat/models_original.py
	cp models.py  [path_to_deeprvat]/deeprvat/deeprvat/deeprvat/
	snakemake -j 60 --snakefile explain.snakefile --profile lsf.profile 
```

Please remember to rename the original `model.py` back once your experiments with SHAP is completed.

```bash
	mv [path_to_deeprvat]/deeprvat/deeprvat/deeprvat/models.py [path_to_deeprvat]/deeprvat/deeprvat/deeprvat/models_shap.py
	mv [path_to_deeprvat]/deeprvat/deeprvat/deeprvat/models_original.py [path_to_deeprvat]/deeprvat/deeprvat/deeprvat/models.py
```