## Mutagenesis experiments

### Generate supplementary figures

You can use [mutagenesis.ipynb](https://github.com/HolEv/deeprvat-analysis/blob/main/feature_importance/mutagenesis/mutagenesis.ipynb) to generate supplementary Figures 3.x and 3.y using [per_pheno_importance.csv](https://github.com/HolEv/deeprvat-analysis/blob/main/feature_importance/mutagenesis/per_pheno_importance.csv) file. 

### Performing mutagenesis experiment

Use `mutate.snakefile` to perform mutagenesis experiments to generate raw files which are used to obtain `per_pheno_importance.csv`. To be able to run these experiments you need to provide the directory of model checkpoints as base_dir in the snakemake file. 

Mutagenesis experiment uses [models.py](https://github.com/bfclarke/deeprvat/blob/master/deeprvat/deeprvat/models.py) file. Make sure to place `explain_mutagenesis_indv.py` script in the same directory as `model.py` and provide the path for `py` in `mutate.snakefile`. 

``
	conda activate [deeprvat_environment]
``
``
	snakemake -j 50 --snakefile mutate.snakefile --profile lsf.profile 
``
