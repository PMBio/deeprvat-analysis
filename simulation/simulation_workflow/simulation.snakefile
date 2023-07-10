from snakemake.utils import Paramspace
from snakemake.utils import min_version
import os

min_version("6.0")

from pathlib import Path
import pandas as pd


configfile: "config.yaml"


debug_flag = config.get("debug", False)
debug = "--debug " if debug_flag else ""

use_sim_causal_variants_only = False
simulation_config_df = pd.read_parquet("simulation_config.parquet")

paramspace = Paramspace(
    simulation_config_df,
    filename_params="*",
    param_sep=":",
)

conda_check = 'conda info | grep "active environment"'

DEEPRVAT_ANALYSIS_DIR = os.environ["DEEPRVAT_ANALYSIS_DIR"]
DEEPRVAT_DIR = os.environ["DEEPRVAT_DIR"]

pipeline_dir = f"{DEEPRVAT_ANALYSIS_DIR}/simulation/simulation_workflow/"
py_sim = f"python {pipeline_dir}"
py_deeprvat = f"python {DEEPRVAT_DIR}/deeprvat"


wildcard_constraints:
    repeat="\d+",


phenotypes = ["sim_phenotype"]

association_testing_maf = config.get("association_testing_maf", 0.01)


n_repeats = config["deeprvat"].get("n_repeats", 6)
n_sim_repeats = 1


##### Define all rule ################################################################
########################################################################################

# if only data sets should be generated (needed for MAF filter experiments (Figure 2b), comment out rule all)

# rule all:
#     input:
#         expand("{params}/eval/metrics.pkl", params=paramspace.instance_patterns),
#         expand(
#             "{params}/eval/variants_per_gene.pdf", params=paramspace.instance_patterns
#         ),
rule all_simulation_datsets:
    input:
        expand(
            "{params}/sim_repeat_{sim_repeat}/simulation_data/phenotype_covariates_simulated.parquet",
            params=paramspace.instance_patterns,
            sim_repeat=range(n_sim_repeats),
        ),
        expand(
            "{params}/eval/variants_per_gene.pdf", params=paramspace.instance_patterns
        ),

##### evauluate simulation ################################################################
########################################################################################


rule evaluate_simulation:
    input:
        expand(
            "{params}/sim_repeat_{sim_repeat}/baseline/{phenotype}/eval/burden_associations.parquet",
            params=paramspace.instance_patterns,
            sim_repeat=range(n_sim_repeats),
            phenotype=phenotypes,
        ),
        expand(
            "{params}/sim_repeat_{sim_repeat}/deeprvat/{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations.parquet",
            params=paramspace.instance_patterns,
            repeat=range(n_repeats),
            sim_repeat=range(n_sim_repeats),
            phenotype=phenotypes,
        ),
    output:
        f"{paramspace.wildcard_pattern}/eval/metrics.pkl",
    params:
        config="config.yaml",
        exp_name=f"{paramspace.wildcard_pattern}",
    threads: 1
    resources:
        mem_mb=16000,
        load=16000,
    shell:
        " && ".join(
        [
            conda_check,
        py_sim + "evaluate_simulation.py evaluate "
                "--n-repeats "
                + str(n_sim_repeats)
                + " {params.exp_name} {params.config}",
            ]
        )


##### run deeprvat ################################################################
########################################################################################

rule all_deeprvat:
    input:
        expand(
            "{params}/sim_repeat_{sim_repeat}/deeprvat/{phenotype}/deeprvat/eval/significant.parquet",
            params=paramspace.instance_patterns,
            phenotype=phenotypes,
            sim_repeat=range(n_sim_repeats),
        ),


module deeprvat_workflow:
    snakefile:
        f"{DEEPRVAT_ANALYSIS_DIR}/phenotype_prediction/cv_deeprvat_training/training_association_testing_with_prefix.snakefile"
    prefix:
        f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/deeprvat"
    config:
        config["deeprvat"]


use rule * from deeprvat_workflow exclude config, choose_training_genes, evaluate, train_bagging, regress, compute_burdens, best_bagging_run, cleanup_burden_cache as deeprvat_*


use rule config from deeprvat_workflow as deeprvat_config with:
    input:
        config=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/{{phenotype}}/deeprvat/config.yaml",
        baseline=expand(
            "{this_pattern}/sim_repeat_{{sim_repeat}}/baseline/{{phenotype}}/eval/burden_associations.parquet",
            this_pattern=paramspace.wildcard_pattern,
        ),
    params:
        baseline_results=lambda wildcards, input: "".join(
            [f"--baseline-results {b} " for b in input.baseline]
        ),


use rule train from deeprvat_workflow as deeprvat_train with:
    params:
        prefix=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/deeprvat",
        phenotypes = " ".join(
            [f"--phenotype {p} "
                f"{p}/deeprvat/input_tensor.zarr "
                f"{p}/deeprvat/covariates.zarr "
                f"{p}/deeprvat/y.zarr"
                for p in ['sim_phenotype']]),


use rule regress from deeprvat_workflow as deeprvat_regress with:
    params:
        prefix=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/deeprvat",


use rule compute_burdens from deeprvat_workflow as deeprvat_compute_burdens with:
    params:
        prefix=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/deeprvat",


use rule best_training_run from deeprvat_workflow as deeprvat_best_training_run with:
    params:
        prefix=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/deeprvat",


use rule cleanup_burden_cache from deeprvat_workflow as deeprvat_cleanup_burden_cache with:
    params:
        prefix=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/deeprvat",



rule all_deeprvat_config:
    input:
        expand(
            "{params}/sim_repeat_{sim_repeat}/{phenotype}/deeprvat/config.yaml",
            sim_repeat=range(n_sim_repeats),
            phenotype=phenotypes,
            params=paramspace.instance_patterns,
        ),



rule config_template_deeprvat:
    input:
        sim_phenotype_file=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/simulation_data/phenotype_covariates_simulated.parquet",
        config="configs/deeprvat_config.yaml",
    output:
        config=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/{{phenotype}}/deeprvat/config.yaml",
    shell:
        " && ".join(
        [
            conda_check,
        py_sim + "simulation_utils.py make-deeprvat-config "
                " --phenotype {wildcards.phenotype}"
                + " --sim-phenotype-file {input.sim_phenotype_file}"
                + " --correction-method FDR --n-training-genes 40"
                + " {input.config} {output.config}",
            ]
        )

###### run baseline #####################################
############################################
rule all_baseline:
    input:
        expand(
            "{params}/sim_repeat_{sim_repeat}/baseline/{phenotype}/eval/burden_associations.parquet",
            params=paramspace.instance_patterns,
            sim_repeat=range(n_sim_repeats),
            phenotype=phenotypes,
        ),

module baseline_workflow:
    snakefile:
        f"{DEEPRVAT_DIR}/pipelines/seed_gene_discovery.snakefile"
    prefix:
        f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/baseline"
    config:
        config["baseline"]


use rule * from baseline_workflow exclude config as baseline_*


rule baseline_config:
    input:
        config="configs/baseline_config.yaml",
        sim_phenotype_file=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/simulation_data/phenotype_covariates_simulated.parquet",
    output:
        f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/baseline/{{phenotype}}/{{vtype}}/config.yaml",
    params:
        rare_maf=str(association_testing_maf),
    threads: 1
    resources:
        mem_mb=1024,
        load=1000,
    shell:
        " && ".join(
        [
            conda_check,
            (
        "seed_gene_pipeline update-config "
                    + "--phenotype {wildcards.phenotype} "
                    + "--variant-type {wildcards.vtype} "
                    + "--rare-maf "
                    + "{params.rare_maf} "
                    + "--simulated-phenotype-file "
                    + "{input.sim_phenotype_file} "
                    + "--maf-column UKB_AF "
                    + " {input.config} "
                    + "{output}"
                ),
            ]
        )


############################################
## spread config

rule spread_config:
    input:
        config="config.yaml",
    output:
        simulation="configs/simulation_config.yaml",  #is bascially not used
        baseline="configs/baseline_config.yaml",
        deeprvat="configs/deeprvat_config.yaml",
    threads: 1
    resources:
        mem_mb=1024,
        load=1000,
    shell:
        " && ".join(
        [
            conda_check,
        py_sim + "simulation_utils.py spread-config "
                "-m baseline -m simulation -m deeprvat "
                " {input.config} configs/",
            ]
        )


# ###### simulate dataset ########################
# ############################################


rule evaluate_simulation_dataset:
    input:
        expand(
            "{params}/sim_repeat_{sim_repeat}/simulation_data/simulated_causal_variants.pkl",
            params=paramspace.instance_patterns,
            sim_repeat=range(n_sim_repeats),
        ),
    output:
        f"{paramspace.wildcard_pattern}/eval/variants_per_gene.pdf",
    resources:
        mem_mb=16000,
        load=16000,
    params:
        config="config.yaml",
        exp_name=f"{paramspace.wildcard_pattern}",
    shell:
        " && ".join(
        [
            conda_check,
        py_sim + "evaluate_simulation.py analyse-variant-properties "
                "--n-repeats "
                + str(n_sim_repeats)
                + " {params.exp_name} {params.config}",
            ]
        )


rule simulate_dataset:
    input:
        config="config.yaml",
        #simulation parameters are read from exp_config param, so the config file can be the same for all 
    output:
        sim_phenotype=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/simulation_data/phenotype_covariates_simulated.parquet",
        sim_causal_variants=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/simulation_data/simulated_causal_variants.pkl",
    threads: 1
    priority: 50
    params:
        exp_config=paramspace.instance,
        debug=debug_flag,
        out_path=f"{paramspace.wildcard_pattern}/sim_repeat_{{sim_repeat}}/simulation_data/",
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * (attempt + 1),
        load=16000,
    script:
        "simulate_data.py"
