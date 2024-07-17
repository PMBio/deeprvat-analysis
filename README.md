# deeprvat-analysis
Code for generating all plots for the paper 'Integration of variant annotations using deep set networks boosts rare variant association genetics' (building on https://github.com/PMBio/deeprvat/)

To run the code, the DeepRVAT conda environment is needed and DeepRVAT has to be installed. Additionaly, the R environment (`r-env`) should be installed from `r-env.yaml`.

The UK Biobank genotype and phenotype data required to run all the analyses presented here is not publicly available. We therefore describe for each experiment the required input data/experiment folder set-up in the corresponding README. For detailed information on the input data, the documentation of the [DeepRVAT](https://github.com/PMBio/deeprvat/) repository should be checked. 

**Instructions to recreate all main figures and related supplementary figures can be found in the respective `figureX` directory.**


The following environment variables need to be exported when using any pipeline in this repository:
- `DEEPRVAT_ANALYSIS_DIR={path_to_this_repository}`
- `DEEPRVAT_DIR={path_to_the_deeprvat_repository}`


