from ast import dump
from email.policy import default
import pandas as pd
import yaml
import os
import sys
from typing import Optional

# import pickle
import logging
import click
from torch.utils.data import DataLoader, Dataset, Subset


logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


@click.group()
def cli():
    pass


@cli.command()
@click.option("--variant-type", type=str)
@click.option("--rare-maf", type=float)
@click.option("--rare-maf-col", type=str, default="UKB_AF")
@click.argument("old_config_file", type=click.Path(exists=True))
@click.argument("new_config_file", type=click.Path())
def update_config(
    old_config_file: str,
    variant_type: Optional[str],
    rare_maf: Optional[float],
    rare_maf_col: Optional[str],
    new_config_file: str,
):
    with open(old_config_file) as f:
        config = yaml.safe_load(f)

    if rare_maf is not None:
        config["data"]["dataset_config"]["min_common_af"][rare_maf_col] = rare_maf
        config["data"]["dataset_config"]["rare_embedding"]["config"]["thresholds"][
            rare_maf_col
        ] = f"{rare_maf_col} < {rare_maf} and {rare_maf_col} > 0"

    with open(new_config_file, "w") as f:
        yaml.dump(config, f)


@cli.command()
@click.option("--module", "-m", multiple=True)
@click.option(
    "--gene-file",
    type=click.Path(exists=True),
    default="genes.parquet",
)
@click.option("--rare-maf-col", type=str, default="UKB_AF")
@click.argument("input_config", type=click.Path(exists=True))
@click.argument("out_path", type=click.Path())
def spread_config(input_config, out_path, module, gene_file, rare_maf_col):
    data_modules = module

    with open(input_config) as f:
        config = yaml.safe_load(f)

    association_maf = config.get("association_testing_maf", 0.01)
    logger.info(
        f"MAF used in association testing for DeepRVAT and baseline: {association_maf} "
    )

    for module in ["baseline", "deeprvat"]:
        data_slots = ["data", "training_data"] if module == "deeprvat" else ["data"]
        for data_slot in data_slots:
            config[module][data_slot]["dataset_config"]["gene_file"] = gene_file
            config[module][data_slot]["dataset_config"]["rare_embedding"]["config"][
                "gene_file"
            ] = gene_file
            # if module == 'baseline':
            #     config[module]['rare_maf'] = association_maf
            if (module == "deeprvat") & (data_slot == "data"):
                # only association testing MAF is changed, training MAF stays at 0.01
                # these columns are updated for the baseline as well but by baseline_scoretest_smk.py update-config
                # based on rare_maf passed in the snakefile
                config[module][data_slot]["dataset_config"]["min_common_af"][
                    rare_maf_col
                ] = association_maf
                config[module][data_slot]["dataset_config"]["rare_embedding"]["config"][
                    "thresholds"
                ][
                    rare_maf_col
                ] = f"{rare_maf_col} < {association_maf} and {rare_maf_col} > 0"

    for data_module in data_modules:
        logger.info(f"Writing config for module {data_module}")
        if data_module not in ["simulation", "baseline", "deeprvat"]:
            ValueError("Invalid data_module")

        with open(f"{out_path}/{data_module}_config.yaml", "w") as f:
            yaml.dump(config[data_module], f)


@cli.command()
@click.option("--phenotype", type=str, default="sim_phenotype")
@click.option("--sim-phenotype-file", type=click.Path(exists=True))
@click.option("--correction-method", type=str, default="FDR")
@click.option("--n-training-genes", type=int, default=30)
@click.argument("input_config", type=click.Path(exists=True))
@click.argument("out_path", type=click.Path())
def make_deeprvat_config(
    input_config,
    out_path,
    phenotype,
    sim_phenotype_file,
    correction_method,
    n_training_genes,
):
    with open(input_config) as f:
        config = yaml.safe_load(f)

    config["phenotypes"] = {
        phenotype: {
            "correction_method": correction_method,
            "n_training_genes": n_training_genes,
        }
    }

    # use phenotype and covariates from simulated phenotype file
    config["data"]["dataset_config"]["sim_phenotype_file"] = sim_phenotype_file
    config["training_data"]["dataset_config"]["sim_phenotype_file"] = sim_phenotype_file

    # make sure simulated  covariates are used (for y_phenotype this is updated automatically based on config['phenotypes'])
    config["data"]["dataset_config"]["x_phenotypes"] = ["sim_cov_1", "sim_cov_2"]
    config["training_data"]["dataset_config"]["x_phenotypes"] = [
        "sim_cov_1",
        "sim_cov_2",
    ]

    # n_repeats of deeprvat model
    if config.get("n_repeats", None) is None:
        config["n_repeats"] = 6
    if config.get("alpha", None) is None:
        config["alpha"] = 0.05

    with open(f"{out_path}", "w") as f:
        yaml.dump(config, f)


if __name__ == "__main__":
    cli()
