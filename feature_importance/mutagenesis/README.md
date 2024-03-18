## Mutagenesis experiments

### Generate supplementary figures

You can use [mutagenesis.ipynb](https://github.com/PMBio/deeprvat-analysis/blob/main/feature_importance/mutagenesis/mutagenesis.ipynb) to generate supplementary Supp. Fig. 3.5 (a) using [per_pheno_importance.csv](https://github.com/PMBio/deeprvat-analysis/blob/main/feature_importance/mutagenesis/per_pheno_importance.csv) file. 

### Perform mutagenesis experiment

Use `mutate.snakefile` to perform mutagenesis experiments to generate raw files which are used to obtain `per_pheno_importance.csv`. 

To be able to run this pipeline, you need to replace the following variables in `mutate.snakefile` file:
* `[path_to_deeprvat]` with the path to your clone of [deeprvat](https://github.com/PMBio/deeprvat/) repository  
* `[path_to_deeprvat-analysis]`  with the path to your clone of this repository  
* `[path_to_inputs]` with the path that contain the outputs of complete model training where for each phenotype following four files are generated:
  * `input_tensor_file = [path_to_inputs]/{pheno}/deeprvat/input_tensor.zarr`
  * `covariates_file = [path_to_inputs]/{pheno}/deeprvat/covariates.zarr`
  * `y_file = [path_to_inputs]/{pheno}/deeprvat/y.zarr`
  * `training_gene_file = [path_to_inputs]/{pheno}/deeprvat/seed_genes.parquet`


#### Try the pipeline on example data

You can try `mutate.snakefile` pipeline on an example data if you already executed the [training and association pipeline on the example data](https://github.com/PMBio/deeprvat/tree/main#try-the-full-training-and-association-testing-pipeline-on-some-example-data) provided in [deeprvat](https://github.com/PMBio/deeprvat/) repository.  


#### Note

Mutagenesis experiment uses the original [models.py](https://github.com/PMBio/deeprvat/blob/master/deeprvat/deeprvat/models.py) file.


```bash
	conda activate [deeprvat_environment]
	cd [path_to_deeprvat-analysis]/feature_importance/mutagenesis
	snakemake -j 60 --snakefile mutate.snakefile --profile lsf.profile 
```
