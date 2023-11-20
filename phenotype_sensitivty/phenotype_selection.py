import sys
import logging
import numpy as np
import math
from pathlib import Path
import click
import yaml

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)

logger = logging.getLogger(__name__)

@click.command()
@click.option("--exp_dir", type=click.Path(exists=True), default=".")
@click.option("--folds", type=int, default= 0)
@click.option("--downsample_percent", type=float, default= 0.1, help="Percentage of phenotypes to remove each fold")
def phenotype_selection(
    exp_dir: str,
    folds: int,
    downsample_percent: float, 
):
    logger.info(f"Downsampling Phenotype Selection Across Folds")
    logger.info(f"   Downsampling WITHOUT Replacement")
    rng = np.random.default_rng(seed=42)

    with open(f'{exp_dir}/base/config.yaml') as f:
        config = yaml.safe_load(f)
    
    selection_list = list(config["phenotypes"].keys())
    num_phenos_to_remove = math.floor(downsample_percent * len(selection_list))
    
    for fold in range(folds):
        
        assert len(selection_list) >= num_phenos_to_remove, \
            "Requested too many folds to sample-without-replacement for the total number of phenotypes."
        
        rng.shuffle(selection_list)
        
        with open(f'{exp_dir}/fold_{fold}/config.yaml') as f:
            foldconfig = yaml.safe_load(f)

        hold_out = {}
        for i in range(num_phenos_to_remove):
            hold_out[selection_list[i]] = foldconfig["phenotypes"].pop(selection_list[i])

        #Sampling WITHOUT replacement.
        selection_list = list({*selection_list} - {*hold_out})

        with open(f'{exp_dir}/fold_{fold}/config.yaml', "w") as f:
            yaml.dump(foldconfig, f)

        with open(f'{exp_dir}/fold_{fold}/held_out_phenos.yaml', "w") as f:
            yaml.dump(hold_out, f)

        logger.info(f"   Fold {fold} complete:")
        logger.info(f"      Held out phenotypes {hold_out.keys()}")


if __name__ == "__main__":
    phenotype_selection()