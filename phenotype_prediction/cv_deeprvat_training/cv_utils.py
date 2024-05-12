import pandas as pd
import yaml
import os
import sys
from typing import Optional
import re
# import pickle
import logging
import click
import copy


logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)
DATA_SLOT_DICT = {
    "deeprvat": ["data", "training_data"],
    "seed_genes": ["data"],
    "alternative_burdens": ["alt_burdens_data"],
}


@click.group()
def cli():
    pass



module_folder_dict = {
    "seed_genes": "baseline",
    "deeprvat": "deeprvat",
    "alternative_burdens": "alternative_burdens",
}


@cli.command()
@click.option("--module", "-m", multiple=True)
@click.option("--fold", type=int)
@click.option("--n-folds", type=int, default=5)
@click.argument("input_config", type=click.Path(exists=True))
@click.argument("out_path", type=click.Path(), default="./")
def spread_config(
    input_config, out_path, module, fold, n_folds
):
    data_modules = module

    with open(input_config) as f:
        config_template = yaml.safe_load(f)
    split = "train"
    cv_path = f"{config_template['cv_path']}/{n_folds}_fold"
    for module in data_modules:
        config = copy.deepcopy(config_template)
        data_slots = DATA_SLOT_DICT[module]
        for data_slot in data_slots:
            sample_file = f"{cv_path}/samples_{split}{fold}.pkl"
            logger.info(f"setting sample file {sample_file}")
            config[data_slot]["dataset_config"]["sample_file"] = sample_file

        if (module == "deeprvat") | (module == "deeprvat_pretrained"):
            logger.info("Writing baseline directories")
            old_baseline = copy.deepcopy(config["baseline_results"])
            logger.info(config["baseline_results"])
        logger.info(f"Writing config for module {module}")
        with open(f"{out_path}/{module_folder_dict[module]}/config.yaml", "w") as f:
            yaml.dump(config, f)




@cli.command()
@click.option("--fold", type=int)
@click.option("--n-folds", type=int, default=5)
@click.argument("input_config", type=click.Path(exists=True))
@click.argument("out_file", type=click.Path())
def generate_test_config(input_config, out_file, fold, n_folds):
    with open(input_config) as f:
        config = yaml.safe_load(f)
    cv_path = f"{config['cv_path']}/{n_folds}_fold"
    split = "test"
    sample_file = f"{cv_path}/samples_{split}{fold}.pkl"
    logger.info(f"setting sample file {sample_file}")
    for data_slot in DATA_SLOT_DICT["deeprvat"]:
        config[data_slot]["dataset_config"]["sample_file"] = sample_file
    with open(out_file, "w") as f:
        yaml.dump(config, f)


if __name__ == "__main__":
    cli()
