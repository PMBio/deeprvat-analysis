## SHAP experiments

### Generate supplementary figures

You can use [explainer.ipynb](https://github.com/PMBio/deeprvat-analysis/blob/main/feature_importance/mutagenesis/mutagenesis.ipynb) to generate Supp  Figure 3.5 (b) and (c) using [shap_relative_per_sampling.pkl](https://github.com/PMBio/deeprvat-analysis/blob/main/feature_importance/shap/shap_relative_per_sampling.pkl)  and [shap_score_per_sampling.pkl](https://github.com/PMBio/deeprvat-analysis/blob/main/feature_importance/shap/shap_score_per_sampling.pkl) files. 

### Running [SHAP DeepExplainer](https://github.com/slundberg/shap)

Use `explain.snakefile` to run SHAP experiments across different samplings to generate raw files which are used to obtain the pickle files above. To be able to run these experiments you need to provide the directory of model checkpoints as base_dir in the snakemake file. 

SHAP experiment uses a modified version [models.py](https://github.com/PMBio/deeprvat/blob/master/deeprvat/deeprvat/models.py) file to accommadate the design of SHAP Deep Explainer. Therefore in order to run SHAP you need to replace the model  with the one provided here models.py.  Make sure to place `explain_shap.py` script in the same directory as `model.py` and provide the path for `py` in `explain.snakefile`. 

```bash
	conda activate [deeprvat_environment]
	snakemake -j 50 --snakefile explain.snakefile --profile lsf.profile 
```
