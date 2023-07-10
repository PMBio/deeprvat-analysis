# Simulation

The simulation pipeline involves the following steps. 

1. Simulating data 
    - For this, we use real genotypes from the UK Biobank and a fixed set of randomly chosen simulated causal genes
    - The simulator then samples causal variants based on the variant annotations, derives gene scores from the causal variants, which are eventually aggregated into phenotype values  (also adding covariates and noise)
    - The pipeline the outputs the simulated phenotypes and covariates   
    - For each simulation repeat, a single phenotype is simulated
2. Running the seed gene discovery/baseline methods on the simulated data
3. Running DeepRVAT on a single phenotype using the seed genes from 2. 
4. Evaluating DeepRVAT and the seed gene discovery methods (power for detectin simulated causal genes)

## Experiment set-up
In totatl, 4 experiment directories for the simulations have to be created. 
Each simulation directory has to provide the same input data as required for the DeepRVAT pipeline. 
In addition, a file mapping all variant ids to their ensembl ids (as exemplarly shown in `../../data/vars_to_ensgid_mapping.parquet`) and a list of genes that should be simulated as causal has to be provided (`../../data/sim_causal_genes.parquet`).

### Supp.Fig 2.3a (including default parameters)
To run the simulations for Supp. Figure 2.3a (which also comprises the default simulation setup used for the calibration experiments), 
the config_2.3a.yaml has to be put into the experiment directory as `config.yaml`, together with `simulation_config-supp_figure_2.3a.parquet` as `simulation_config.parquet`.
Then, the pipeline can be run using `simulation.snakefile`

### Figure 2b
For this experiment, we simulated datasets as defined in `simulation_config-figure_2b.parquet`. Then, the seed gene pipeline and DeepRVAT are applied to each simulated data set three times, every time using a different threshold for the variant maf (1%, 0.1%, 0.01%).
Technically, this means that 3 experiments have to be crated, where  `config.yaml` files only vary in the  `association_testing_maf` value. 
To only simulate the data set in one of the 3 directories, comment out `rule all` and comment in `rule all_simulated_data`. 
Once the data set has been simulated in one of the 3 directories, it has to be copied to the remaining ones and the association testing has to be run in all 3 experiments using `rule all`. 


## Analyzing the simulations to generate the paper figures

To analyze the simulation experiments, `simulation_analysis/analyse_simulations.py` has to be run. This script uses `simulation_analysis/config_eval.yaml`, where the paths to the simulation experiment directories have to be set accordingly. 
Once the script is done, the paper figures can be generated using `paper_figures/fig_2_simulation_plots.Rmd`
