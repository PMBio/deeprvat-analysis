import sys
import logging
import numpy as np
import math
from pathlib import Path
import click
import yaml
import pandas as pd

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)

logger = logging.getLogger(__name__)

@click.command()
@click.option("--exp_dir", type=click.Path(exists=True), default=".")
@click.option("--reps", type=int, default= 0)
@click.option("--downsample_percent", type=float, default= 0.1, help="Percentage of phenotypes to remove in each rep")
def phenotype_selection(
    exp_dir: str,
    reps: int,
    downsample_percent: float, 
):
    '''
    This script is used for running the phenotype-downsample analysis. 
    Works with the cv_split model setup.
    '''
    logger.info(f"Downsampling Phenotype Selection Across reps")
    logger.info(f"   Downsampling WITHOUT Replacement")
    rng = np.random.default_rng(seed=42)

    with open(f'{exp_dir}/base/deeprvat_config.yaml') as f:
        config = yaml.safe_load(f)
    
    phenotypes = config["training"]["phenotypes"]
    all_phenos = config["phenotypes"]
    num_phenos_to_remove = math.floor(downsample_percent * len(phenotypes))
    association_only_phenos = set(all_phenos)-set(phenotypes)
    logger.info(f"  Phenotypes already reserved for association testing: {association_only_phenos}")

    pheno_seed_genes = {}
    for pheno in phenotypes:
        seed_gene_df = pd.read_parquet(f'{exp_dir}/base/cv_split0/deeprvat/{pheno}/deeprvat/seed_genes.parquet', engine="pyarrow")
        pheno_seed_genes[pheno] = len(seed_gene_df)

        del seed_gene_df

    #rank phenotypes in descending order of seed genes
    ranked_phenos_descending = sorted(pheno_seed_genes.items(), key=lambda kv:kv[1], reverse=True)

    groups = {}
    prev = 0
    step = len(phenotypes) // num_phenos_to_remove
    for i,idx in enumerate(range(step,len(phenotypes),step)):
        groups[f"group_{i}"]=ranked_phenos_descending[prev:idx]
        prev = idx
    #re-assign phenotypes that aren't evenly divided in previously. Assign begining with the last group.
    if len(phenotypes) % step != 0:
        mod = len(phenotypes)%step
        #reassign phenotypes, beginning with the last group
        g = len(groups) -1  #subtract 1 for 0 indexing
        for i in range(1,mod+1):
            groups[f"group_{g}"].append(ranked_phenos_descending[-i])
            g -= 1 
    logger.info(f" Groups:\n{groups}")

    assert len(groups) == num_phenos_to_remove
    assert len(groups['group_0']) >= reps, \
            "Requested too many reps to sample-without-replacement for the total number of phenotypes."
    
    hold_out = {f"rep_{rep}" : [] for rep in range(reps)}
    for rep in range(reps):
        rep_select = groups
        
        with open(f'{exp_dir}/base/deeprvat_config.yaml','r') as f:
            repconfig = yaml.safe_load(f)

        for i in range(num_phenos_to_remove):
            remove_pheno = rep_select[f'group_{i}'][0][0] #always just select idx 0, as we shuffle each loop

            repconfig["training"]["phenotypes"].remove(remove_pheno)
            hold_out[f"rep_{rep}"].append(remove_pheno)

            #Sample WITHOUT replacement
            rep_select[f'group_{i}'].pop(0)
            rng.shuffle(rep_select[f'group_{i}'])

        with open(f'{exp_dir}/rep_{rep}/deeprvat_config.yaml', "w") as f:
            yaml.dump(repconfig, f) #save new set of training_phenotypes

        logger.info(f"   rep {rep} complete:")
        logger.info(f"      Held out phenotypes {hold_out[f'rep_{rep}']}")

    with open(f'{exp_dir}/held_out_phenos.yaml', "w") as f:
        yaml.dump(hold_out, f)


if __name__ == "__main__":
    phenotype_selection()