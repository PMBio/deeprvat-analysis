import sys
import logging
from typing import Optional
import numpy as np
import pandas as pd
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
@click.option("--seed_gene_file", type=str, default="seed_genes.parquet")
@click.option("--folds", type=int, default= 0, help="The number of different seed-gene setups, i.e. number of experiments each with a different set of seed genes.")
@click.option("--downsample_percent", type=float, default= 0.1, help="Percentage of total seed genes to remove.")
@click.option("--min_keep_percent", type=float, default= 0.5, help="At least this percentage of seed genes per phenotype must be kept.")
@click.option("--min_seed_genes", type=int, default= 4)
@click.option("--cv_splits", type=int, default= 5)
def seed_gene_selection(
    exp_dir: str,
    seed_gene_file: Optional[str],
    folds: int,
    downsample_percent: float, 
    min_keep_percent: float, 
    min_seed_genes: int,
    cv_splits: 5,
):
    '''
    This script is used for running the seed-gene-downsample analysis. 
    Works with the cv_split model setup.
    '''
    logger.info(f"Downsampling Seed Gene Selection")
    rng = np.random.default_rng(seed=42)

    assert (1-min_keep_percent) > downsample_percent
    
    sg_dict = {}
    total_seed_genes = 0
    resample_phenos = []

    with open(f'{exp_dir}/base/deeprvat_config.yaml') as f:
        config = yaml.safe_load(f)
    
    phenotypes = config["training"]["phenotypes"]
    all_phenos = config["phenotypes"]
    association_only_phenos = set(all_phenos)-set(phenotypes)
    logger.info(f"  Phenotypes reserved for association testing: {association_only_phenos}")
    
    for pheno in phenotypes:
        seed_gene_df = pd.read_parquet(f'{exp_dir}/base/cv_split0/deeprvat/{pheno}/deeprvat/{seed_gene_file}', engine="pyarrow") ##all cv_splits have the same seed_gene.parquet for a given phenotype
        sg_dict[pheno] = len(seed_gene_df)
        if sg_dict[pheno] <= min_seed_genes:
            logger.info(f"  Too few seed genes. Keeping all seed genes for {pheno}. Number of seed genes = {sg_dict[pheno]}")
        else: 
            resample_phenos.append(pheno)
        total_seed_genes += len(seed_gene_df)
        del seed_gene_df
    total_genes_remove = math.ceil(downsample_percent * total_seed_genes)
    logger.info(f"Total Genes to Remove = {total_genes_remove}")
    classes = len(resample_phenos)

    # Solutions will hover closely to evenly split selections for each phenotype
    remove_quantity = rng.multinomial(total_genes_remove, [1/classes]*classes, size=folds) # shape = (folds, resample_phenos)
    sg_remove_dict = {pheno: remove_quantity[:,col] for col, pheno in enumerate(resample_phenos)}

    reassign = {fold: 0 for fold in range(folds)}
    locked_phenos = {fold: [] for fold in range(folds)} 
    for fold in range(folds):
        for pheno in resample_phenos: 
            up_bound = math.floor((1-min_keep_percent)*sg_dict[pheno])
            min_count = sg_dict[pheno] - min_seed_genes 

            if sg_remove_dict[pheno][fold] > min(up_bound, min_count):
                while sg_remove_dict[pheno][fold] > min(up_bound, min_count):  
                    logger.info(f"  Check 1 : Selected too many genes to remove from {pheno}. Reducing by 1 and trying again.")
                    sg_remove_dict[pheno][fold] -= 1
                    reassign[fold] += 1

                locked_phenos[fold].append(pheno)
    
    #Randomly reassign the errored knockouts to other phenotypes
    for fold in range(folds):
        rng.shuffle(resample_phenos)
        while reassign[fold] > 0:
            for pheno in resample_phenos:
                if pheno not in locked_phenos[fold]:
                    
                    up_bound = math.floor((1-min_keep_percent)*sg_dict[pheno])
                    min_count = sg_dict[pheno] - min_seed_genes 
                    
                    if sg_remove_dict[pheno][fold] < min(up_bound, min_count): 
                        print(f"  Check 2: Reassigning genes to {pheno} - {fold}.")
                        sg_remove_dict[pheno][fold] += 1
                        reassign[fold] -= 1
                        if reassign[fold] == 0:
                            break
                    
                    if sg_remove_dict[pheno][fold] == min(up_bound, min_count) :
                        locked_phenos[fold].append(pheno)

            if (({*resample_phenos} & {*locked_phenos[fold]}) == {*resample_phenos} ) and (reassign[fold] > 0):  
                #if all phenos are now locked and there are still genes to knockout
                logger.info(f"   Too high kept amount [{min_keep_percent*100}%] in combination with too high downsample amount [{downsample_percent*100}%]")
                logger.info(f"   Could not distribute all portion of genes to be knocked out.")
                logger.info(f"   Try again with lowered min_keep_percent and/or lowered downsample_percent.")
                assert reassign[fold] == 0

    #Select seed genes for each cv_split and save as new parquet file
    for fold in range(folds):
        for split in range(cv_splits):
            for pheno in phenotypes:
                seed_gene_df = pd.read_parquet(f'{exp_dir}/base/cv_split{split}/deeprvat/{pheno}/deeprvat/{seed_gene_file}', engine="pyarrow")

                with open(f'{exp_dir}/base/cv_split{split}/deeprvat/{pheno}/deeprvat/config.yaml') as f:
                    baseconfig = yaml.safe_load(f)

                p = Path(f'{exp_dir}/sg_set_{fold}/cv_split{split}/deeprvat/{pheno}/deeprvat')
                p.mkdir(parents=True,exist_ok=True)

                if (pheno not in resample_phenos) or (sg_remove_dict[pheno][fold] == 0):
                    seed_gene_df.to_parquet(f"{exp_dir}/sg_set_{fold}/cv_split{split}/deeprvat/{pheno}/deeprvat/seed_genes.parquet", engine="pyarrow")
                    logger.info(f"    sg_set{fold} cv_split{split} - {pheno}: No seed genes removed. Number seed genes = {len(seed_gene_df)}")
                else:
                    logger.info(f"    sg_set{fold} cv_split{split} - {pheno}: Number of seed genes BEFORE selection = {len(seed_gene_df)}")
                    nrows = len(seed_gene_df) - sg_remove_dict[pheno][fold] #number of seed genes to keep
                    seed_gene_new_df = seed_gene_df.sample(nrows) #keep n-number of rows randomly
                    logger.info(f"     Removing {sg_remove_dict[pheno][fold]} seed genes from {pheno} for cv_split {split} in sg_set_{fold}.")
                    logger.info(f"     Resulting number of seed genes = {len(seed_gene_new_df)}")
                    seed_gene_new_df.to_parquet(f"./sg_set_{fold}/cv_split{split}/deeprvat/{pheno}/deeprvat/seed_genes.parquet", engine="pyarrow")
                    logger.info(f"     sg_set{fold} cv_split{split} - {pheno}: Saved seed_genes.parquet")
                    del seed_gene_new_df
                
                with open(f'{exp_dir}/sg_set_{fold}/cv_split{split}/deeprvat/{pheno}/deeprvat/config.yaml', "w") as f:
                    yaml.dump(baseconfig, f)
                    
                del seed_gene_df

            if association_only_phenos:
                for association_pheno in association_only_phenos:
                    with open(f'{exp_dir}/base/cv_split{split}/deeprvat/{association_pheno}/deeprvat/config.yaml') as f:
                        baseconfig = yaml.safe_load(f)

                    p = Path(f'{exp_dir}/sg_set_{fold}/cv_split{split}/deeprvat/{association_pheno}/deeprvat')
                    p.mkdir(parents=True,exist_ok=True)
                    with open(f'{exp_dir}/sg_set_{fold}/cv_split{split}/deeprvat/{association_pheno}/deeprvat/config.yaml', "w") as f:
                        yaml.dump(baseconfig, f)


if __name__ == "__main__":
    seed_gene_selection()