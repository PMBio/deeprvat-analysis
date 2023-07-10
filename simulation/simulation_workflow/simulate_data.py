from ast import dump
from email.policy import default
import pandas as pd
import yaml
import os
import sys

# import pickle
from simulator import SimulatedPhenotypes
import logging
import click
from torch.utils.data import DataLoader, Dataset, Subset


logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

column_renamer_reverse = {
    "vprop": "prop_causal_variants",
    "maxaf": "max_rare_af",
    "vnoise": "var_noise",
    "vgen": "var_genetic_effects",
    "vvarnoise": "var_expl_variant_noise",
    "vcont": "var_expl_cont",
    "vbin": "var_expl_binary",
    "vmaf": "var_expl_maf",
    "vextra": "var_extra_group",
    "expbin": "binary_anno_weight_exp",
}


def simulate_data_(
    config,
    debug,
    input_tensor_file,
    prop_causal_variants: float,
    max_rare_af: float,
    var_expl_variant_noise: float,
    var_expl_binary: float,
    var_expl_maf: float,
    var_expl_cont: float,
    var_genetic_effects: float,
    var_noise: float,
    var_extra_group: float = None,
    binary_anno_weight_exp: int = 1,
):
    logger.info("simulating dataset")

    simulator = SimulatedPhenotypes(
        input_tensor_file=input_tensor_file,
        debug=debug,
        config=config,
        prop_causal_variants=prop_causal_variants,
        max_rare_af=max_rare_af,
        var_expl_variant_noise=var_expl_variant_noise,
        var_expl_maf=var_expl_maf,
        var_expl_binary=var_expl_binary,
        var_expl_cont=var_expl_cont,
        var_genetic_effects=var_genetic_effects,
        var_noise=var_noise,
        var_extra_group=var_extra_group,
        binary_anno_weight_exp=binary_anno_weight_exp,
        **config["simulation_config"]["config"],
    )

    simulator.simulate_data_set()

    sim_pheno_and_covariates = simulator.phenotype_df
    sim_causal_genes = simulator.causal_genes
    sim_causal_variants = simulator.all_causal_var_ids
    variant_scores = simulator.var_scores

    return (
        sim_pheno_and_covariates,
        sim_causal_genes,
        sim_causal_variants,
        variant_scores,
    )


def simulate_data(debug: bool, config_file: str, simulation_params, out_dir: str):
    if "vmaf" in simulation_params.keys():
        logger.info(f"Renaming short simulation params: {simulation_params}")
        new_params = {}
        for key in simulation_params.keys():
            this_val = simulation_params[key]
            new_params[column_renamer_reverse[key]] = this_val
        simulation_params = new_params
        logger.info(f"Simulation params after renaming: {simulation_params}")

    with open(config_file) as f:
        config = yaml.safe_load(f)["simulation"]
    input_tensor_file = f"{out_dir}/input_tensor.zarr"
    (
        sim_pheno_and_covariates,
        sim_causal_genes,
        sim_causal_variants,
        variant_scores,
    ) = simulate_data_(config, debug, input_tensor_file, **simulation_params)

    logger.info("Writing output files")
    sim_pheno_and_covariates.to_parquet(
        f"{out_dir}/phenotype_covariates_simulated.parquet"
    )
    sim_causal_genes.to_parquet(f"{out_dir}/causal_genes_simulated.parquet")
    variant_scores = (
        variant_scores.reset_index()
        .rename(columns={0: "variant_score"})
        .assign(
            sim_causal=lambda df: df["id"].map(
                lambda x: True if x in sim_causal_variants else False
            )
        )
    )
    variant_scores.to_parquet(f"{out_dir}/variant_scores.parquet")
    with open(f"{out_dir}/simulated_causal_variants.pkl", "wb") as f:
        pickle.dump(sim_causal_variants, f)


if __name__ == "__main__":
    simulate_data(
        snakemake.params["debug"],
        snakemake.input["config"],
        snakemake.params["exp_config"],
        snakemake.params["out_path"],
    )
