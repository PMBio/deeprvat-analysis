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
@click.option("--correction-method", type=str, default="FDR")
@click.option("--n-training-genes", type=int, default=40)
@click.argument("input_config", type=click.Path(exists=True))
@click.argument("out_path", type=click.Path(), default="./")
def spread_config(
    input_config,
    out_path,
    module,
    fold,
    correction_method,
    n_training_genes,
):
    data_modules = module

    with open(input_config) as f:
        config_template = yaml.safe_load(f)

    association_maf = config_template.get("association_testing_maf", 0.001)
    logger.info(
        f"MAF used in association testing for DeepRVAT and baseline: {association_maf} "
    )
    for split in ["train"]:
        for module in data_modules:
            config = copy.deepcopy(config_template)
            # data_slots = ['data', 'training_data'] if module == 'deeprvat' else ['data']
            data_slots = DATA_SLOT_DICT[module]
            for data_slot in data_slots:
                annotations = config[module][data_slot]["dataset_config"]["annotations"]
                af_pattern = re.compile(r".*(_MAF|_AF)\b")
                rare_maf_col = [s for s in annotations if af_pattern.match(s)]
                print(rare_maf_col)
                assert len(rare_maf_col) == 1
                rare_maf_col = rare_maf_col[0]
                print(rare_maf_col)
                config[module][data_slot][
                    "gt_file"
                ] = f"cv_data/genotypes_{split}{fold}.h5"
                config[module][data_slot]["dataset_config"][
                    "phenotype_file"
                ] = f"cv_data/genotypes_{split}{fold}_phenotypes.parquet"
                # if module == 'baseline':
                #     config[module]['rare_maf'] = association_maf
                if (module == "deeprvat") | (module == "deeprvat_pretrained"):
                    logger.info("writing phenotype config")
                    config[module]["phenotypes"] = {
                        phenotype: {
                            "correction_method": correction_method,
                            "n_training_genes": n_training_genes,
                        }
                        for phenotype in config["phenotypes"]
                    }
                    if data_slot == "data":
                        config[module][data_slot]["dataset_config"]["min_common_af"][
                            rare_maf_col
                        ] = association_maf
                        config[module][data_slot]["dataset_config"]["rare_embedding"][
                            "config"
                        ]["thresholds"][
                            rare_maf_col
                        ] = f"{rare_maf_col} < {association_maf} and {rare_maf_col} > 0"
                    logger.info("Writing baseline directories")
                    baseline_base_path = f"{os.getcwd()}/cv_split{fold}/baseline"
                    logger.info(f"Setting baseline path to {baseline_base_path}")
                    config[module]["baseline_results"] = [
                        {"base": baseline_base_path, "type": test_name}
                        for test_name in [
                            "plof/burden",
                            "plof/skat",
                            "missense/burden",
                            "missense/skat",
                        ]
                    ]
                    logger.info(config[module]["baseline_results"])
            logger.info(f"Writing config for module {module}")
            split_suffix = "_test" if split == "test" else ""
            with open(
                f"{out_path}/{module_folder_dict[module]}/config{split_suffix}.yaml",
                "w",
            ) as f:
                yaml.dump(config[module], f)


@cli.command()
@click.option("-m", "--module", type=str, default="deeprvat")
@click.option("--fold", type=int)
@click.argument("input_config", type=click.Path(exists=True))
@click.argument("out_file", type=click.Path())
def generate_test_config(input_config, out_file, module, fold):
    with open(input_config) as f:
        config = yaml.safe_load(f)
    split = "test"
    for data_slot in DATA_SLOT_DICT[module]:
        config[data_slot]["gt_file"] = f"cv_data/genotypes_{split}{fold}.h5"
        config[data_slot]["dataset_config"][
            "phenotype_file"
        ] = f"cv_data/genotypes_{split}{fold}_phenotypes.parquet"
    with open(out_file, "w") as f:
        yaml.dump(config, f)


if __name__ == "__main__":
    cli()
