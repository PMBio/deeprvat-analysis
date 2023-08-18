# from prs_comparison import get_burdens
import pickle
import pandas as pd
import sys
import os
import logging
from pathlib import Path
import click
import yaml
import zarr
import numpy as np
from prs_comparison import PHENOTYPE_MAPPING
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def get_burdens(burden_dir, query_genes, burden_type):
    all_burdens = zarr.open(Path(burden_dir) / "burdens.zarr")  # [:, :, repeat]
    genes = pd.Series(np.load(Path(burden_dir) / "genes.npy"))
    logger.info(f"All burdens shape: {all_burdens.shape}")

    n_missing_genes = len(query_genes) - len(set(genes).intersection(query_genes))
    query_genes = set(genes).intersection(query_genes)

    logger.info(
        f"Number of genes query genes (genebass hits) that "
        f"are not in DeepRVAT results and therefore removed: {n_missing_genes}"
    )

    gene_id_index_dict = {gene: genes[genes == gene].index[0] for gene in query_genes}
    if burden_type == "deeprvat":
        n_repeats = all_burdens.shape[2]
        out_burdens = {
            gene: {repeat: {} for repeat in range(n_repeats)}
            for gene in gene_id_index_dict.keys()
        }
        for gene, index in gene_id_index_dict.items():
            logger.info(f"Retreiving burdens for gene {gene}")
            for repeat in range(n_repeats):
                burdens_this_gene = all_burdens.oindex[:, index, repeat]

                out_burdens[gene][repeat] = burdens_this_gene
    else:
        assert len(all_burdens.shape) == 2
        out_burdens = {}
        for gene, index in gene_id_index_dict.items():
            logger.info(f"Retreiving burdens for gene {gene}")
            burdens_this_gene = all_burdens.oindex[:, index]
            out_burdens[gene] = burdens_this_gene

    return out_burdens


@click.group()
def cli():
    pass


@cli.command()
@click.option("--discoveries-file", type=click.Path(exists=True))
@click.option("--burdens-out", type=click.Path(), default=None)
@click.option("--burden-input-dir", type=click.Path(exists=True))
@click.option("--burden-type", "-s")
@click.option("--debug", is_flag=True)
def extract_burdens(
    burden_input_dir, burdens_out, burden_type, discoveries_file, debug
):
    logger.info(f"Reading discoveries file {discoveries_file}")
    discoveries = pd.read_parquet(discoveries_file)  # TODO change this to 0.4

    logger.info(f"Retrieving for all phenotypes {burden_input_dir}")
    logger.info(f'Considered phenotypes are {discoveries["Trait"].unique()}')

    logger.info("Extracting genes based on FDR threshold")

    query_genes = set(discoveries["gene"])
    logger.info(f"Number of query genes: {len(query_genes)}")

    logger.info(f"Total number of genes: {len(query_genes)}")

    if debug:
        query_genes = list(query_genes)[:2]

    logger.info("Retreiving burdens")
    burden_dict = get_burdens(burden_input_dir, query_genes, burden_type)

    with open(burdens_out, "wb") as f:
        pickle.dump(burden_dict, f)

    logger.info("succesfully completed")


@cli.command()
@click.option("--sig-file", type=click.Path(exists=True), multiple=True)
@click.option("--train-pheno", multiple=True)
@click.argument("out-file", type=click.Path())
def combine_significant(sig_file, train_pheno, out_file):
    all_discoveries = pd.concat([pd.read_parquet(file) for file in sig_file])
    discovery_list = []
    for pheno in all_discoveries['phenotype'].unique():
        if pheno in train_pheno:
            logger.info(f'Using DeepRVAT with seed genes as DeepRVAT for pheno {pheno}')
            this_deeprvat = all_discoveries.query("phenotype == @pheno & Method == 'DeepRVAT'")
        else:
            logger.info(f'Using DeepRVAT without seed genes as DeepRVAT for pheno {pheno}')
            this_deeprvat = all_discoveries.query("phenotype == @pheno & Method == 'DeepRVAT wo baseline'")
            this_deeprvat['Method'] = 'DeepRVAT'
        discovery_list.append(this_deeprvat)
    baseline_discoveries = all_discoveries.query("Method not in ['DeepRVAT', 'DeepRVAT wo baseline']")
    discovery_list.append(baseline_discoveries)
    all_discoveries = pd.concat(discovery_list).rename(columns={'phenotype':'Trait'})
    this_pheno_mapping = dict(zip([i.replace('_', ' ') for i in PHENOTYPE_MAPPING.keys()], 
        PHENOTYPE_MAPPING.values()))
    all_discoveries['Trait'] =  all_discoveries["Trait"].replace(this_pheno_mapping)
    logger.info(f"Number of discoveries:  {all_discoveries.groupby(['Trait', 'Method']).size()}")

    all_discoveries.to_parquet(out_file)


if __name__ == "__main__":
    cli()
