import logging
import os
import random
import sys
import pickle
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import click
import dask.dataframe as dd
import h5py
import math
import numpy as np
import pandas as pd
import yaml
from deeprvat.data import DenseGTDataset
from scipy.sparse import coo_matrix, spmatrix
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


class GotNone(Exception):
    pass



def extract_gene(
    G_full: spmatrix,
    gene: int,
    grouped_annotations: pd.DataFrame,
    annotation_cols: Iterable[str],
) -> Tuple[np.ndarray, np.ndarray]:
    # Find variants present in gene
    # Convert sparse genotype to CSC
    # Slice genotype and annotation matrices by variants present in gene

    annotation_df = grouped_annotations.get_group(gene)
    variant_ids = annotation_df.index.unique().to_numpy()
    assert len(variant_ids) == len(annotation_df)
    # assert set(annotation_df.columns) == annotation_cols
    annotations = annotation_df[annotation_cols].to_numpy()

    ## sparse Compressed Sparse Column
    genotypes = G_full.tocsc()[:, variant_ids].todense()
    # Important
    # Important: cast G into numpy array. Otherwise it will be a matrix and
    # the * operator does matrix mutiplication (.dot()) instead of scalar multiplication (.multiply())
    genotypes = np.asarray(genotypes)

    return genotypes, annotations


def build_data_(
    gene_ids,
    G_full,
    grouped_annotations: pd.DataFrame,
    annotation_cols: Iterable[str],
    genotype_file: h5py.File,
    annotation_file: h5py.File,
):
    for gene in tqdm(gene_ids, file=sys.stdout):
        genotypes, annotations = extract_gene(
            G_full, gene, grouped_annotations, annotation_cols
        )
        genotype_file.create_dataset(str(gene), data=genotypes)
        annotation_file.create_dataset(str(gene), data=annotations)


@click.group()
def cli():
    pass


@cli.command()
@click.option("--phenotype", type=str)
@click.option("--variant-type", type=str)
@click.option("--simulated-phenotype_file", type=str)
@click.argument("old_config_file", type=click.Path(exists=True))
@click.argument("new_config_file", type=click.Path())
def update_config(
    old_config_file: str,
    phenotype: Optional[str],
    simulated_phenotype_file: str,
    variant_type: Optional[str],
    new_config_file: str,
):
    logger.info("Reading existing config file")
    with open(old_config_file) as f:
        config = yaml.safe_load(f)
    if simulated_phenotype_file is not None:
        logger.info("Using simulated phenotype file")
        config["data"]["dataset_config"][
            "sim_phenotype_file"
        ] = simulated_phenotype_file

    if phenotype is not None:
        logger.info("Setting y_phenotype in dataset")
        config["data"]["dataset_config"]["y_phenotypes"] = [phenotype]
        # config['training_data']['dataset_config']['y_phenotypes'] = [phenotype]
    rare_threshold_config = config["data"]["dataset_config"]["rare_embedding"][
        "config"
    ]["thresholds"]

    variant_filter_dict = {
        "missense": "Consequence_missense_variant",
        "plof": "is_plof",
        "disruptive_missense": "disruptive_missense",
        "plof_disruptive_missense": "plof_or_disruptive_missense",
        "synonymous": "Consequence_synonymous_variant",
    }
    if variant_type is not None:
        logger.info(f"Variant type is {variant_type}")
        filter_column = variant_filter_dict[variant_type]
        rare_threshold_config[filter_column] = f"{filter_column} == 1"
    else:
        logger.info(
            "No variant masking specified. Using all variants meeting MAF threshold"
        )

    print(f"Rare threshold config: {rare_threshold_config}")

    logger.info("Writing modified config file")
    with open(new_config_file, "w") as f:
        yaml.dump(config, f)

    logger.info("Done")


def make_dataset_(
    config: Dict,
    pickled_dataset_file: str = None,
    debug: bool = False,
    data_key="data",
) -> Dataset:
    data_config = config[data_key]

    if pickled_dataset_file is not None and os.path.isfile(pickled_dataset_file):
        logger.info("Loading pickled dataset")
        with open(pickled_dataset_file, "rb") as f:
            dataset = pickle.load(f)
    else:
        logger.info("Instantiating dataset")
        logger.info(
            data_config["dataset_config"]["rare_embedding"]["config"]["thresholds"]
        )
        dataset = DenseGTDataset(
            gt_file=data_config["gt_file"],
            skip_y_na=True,
            skip_x_na=True,
            **data_config["dataset_config"],
        )
        logger.info("Writing pickled data set")
        with open(pickled_dataset_file, "wb") as handle:
            pickle.dump(dataset, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if debug:
        logger.info("Debug mode: Using only 1000 samples")
        batch_size = 1000
    else:
        batch_size = len(dataset)

    logger.info(f"read dataset, batch size {batch_size}")
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,  # reading from dataloader config
        collate_fn=dataset.collate_fn,
        **data_config["dataloader_config"],
    )

    logger.info("Loading data")
    data_full = next(iter(dataloader))

    logger.info("Data succesfully loaded. Data set generation completed.")

    return dataset, data_full


@cli.command()
@click.option("--debug", is_flag=True)
@click.option("--data-key", type=str, default="data")
@click.option("--pickled-dataset-file", type=str, default=None)
@click.argument("config-file", type=click.Path(exists=True))
@click.argument("out-file", type=click.Path())
def make_dataset(
    debug: bool,
    data_key: str,
    config_file: str,
    pickled_dataset_file: str,
    out_file: str,
):
    with open(config_file) as f:
        config = yaml.safe_load(f)

    _, ds = make_dataset_(
        config,
        debug=debug,
        data_key=data_key,
        pickled_dataset_file=pickled_dataset_file,
    )

    with open(out_file, "wb") as f:
        pickle.dump(ds, f)


@cli.command()
@click.option("--debug", is_flag=True)
@click.option("--n-chunks", type=int)
@click.option("--chunk", type=int)
@click.option("--sample-file", type=click.Path(exists=True))
@click.option("--dataset-file", type=click.Path(exists=True))
@click.option("--data-file", type=click.Path(exists=True))  # dataset_full
@click.argument("config-file", type=click.Path(exists=True))
@click.argument("genotype_file", type=click.Path())
@click.argument("annotation_file", type=click.Path())
@click.argument("x_file", type=click.Path())
@click.argument("y_file", type=click.Path())
def build_data(
    debug: bool,
    dataset_file: Optional[str],
    data_file: Optional[str],
    sample_file: Optional[str],
    config_file: str,
    genotype_file: str,
    annotation_file: str,
    x_file: str,
    y_file: str,
    n_chunks: Optional[int] = None,
    chunk: Optional[int] = None,
):
    # res_dir = Path(f'{out_dir}/results')
    # Read config
    with open(config_file) as f:
        config = yaml.safe_load(f)

    if data_file is not None and dataset_file is not None:
        logger.info("Loading pickled data")
        with open(dataset_file, "rb") as f:
            dataset = pickle.load(f)
        with open(data_file, "rb") as f:
            data_full = pickle.load(f)
    else:
        logger.info("One of data_file or dataset_file is None; generating data")
        dataset, data_full = make_dataset_(config, debug=debug)
    all_samples = np.array([int(i) for i in data_full["sample"]])
    if sample_file is not None:
        logger.info(f"Using sample file {sample_file}")
        with open(sample_file, "rb") as f:
            samples = pickle.load(f)["training_samples"]
        training_dataset_file = (
            f"{'/'.join(sample_file.split('/')[:-2])}/training_dataset.pkl"
        )
        # load this file to remap the sample ids
        with open(training_dataset_file, "rb") as f:
            ref_training_datset = pickle.load(f)

        # Gett actual sample ids (from dataset the splitting was based on)
        this_sample_ids = ref_training_datset.samples[samples].astype("int").tolist()
        if len(set(this_sample_ids) - set(all_samples)) > 0:
            logger.info(
                "Not all required sample ids are part of the data set \
                Only selecting those samples that overlap"
            )

            this_sample_ids = list(set(this_sample_ids).intersection(set(all_samples)))

        logger.info("Remapping sample indices")

        # Remap to get array ids in our data set

        this_data_idx = [
            np.where(all_samples == this_id)[0][0] for this_id in this_sample_ids
        ]

    else:
        this_data_idx = [i for i in range(len(data_full["sample"]))]

    G_full = data_full["rare_variant_annotations"]
    all_variants = np.unique(G_full.col)  # SparseGenotype
    # all_variants = np.unique(G_full)       #PaddedAnnotations
    if isinstance(G_full, spmatrix):
        G_full = G_full.tocsr()
    G_full = G_full[this_data_idx]
    ## HAKIME
    logger.info(f"G_full shape: {G_full.shape}")
    logger.info(f"all_variants shape: {all_variants.shape}")

    X = data_full["x_phenotypes"].numpy()[this_data_idx]
    X = np.hstack((np.ones((X.shape[0], 1)), X))  # add bias column
    Y = data_full["y"].numpy()[this_data_idx]
    logger.info(f"X shape: {X.shape}")
    logger.info(f"Y shape: {Y.shape}")

    logger.info("Subsetting data into discovery/validation")
    train_proportion = config.get("train_proportion", 1.0)
    logger.info(f"'train' proportion: {train_proportion}")

    logger.info("Grouping variants by gene")
    annotation_cols = config["data"]["dataset_config"]["rare_embedding"]["config"][
        "annotations"
    ]
    annotation_df = dataset.annotation_df.query("id in @all_variants")
    annotation_df = annotation_df[["gene_ids"] + annotation_cols]
    exploded_annotations = (
        dataset.annotation_df.query("id in @all_variants")
        .explode("gene_ids")
        .drop_duplicates()
    )  # row can be duplicated if a variant is assigned to a gene multiple times
    grouped_annotations = exploded_annotations.groupby("gene_ids")
    gene_ids = pd.read_parquet(dataset.gene_file, columns=["id"])["id"].to_list()
    gene_ids = list(
        set(gene_ids).intersection(set(exploded_annotations["gene_ids"].unique()))
    )
    logger.info(f"Number of genes to test: {len(gene_ids)}")

    if debug:
        logger.info("Debug mode: Using only 5000 genes")
        gene_ids = gene_ids[:5000]

    n_total_genes = len(gene_ids)
    logger.info("Reading variant file")
    logger.info("Training split: Running tests for each gene")

    if chunk is not None:
        if n_chunks is None:
            raise ValueError("n_chunks must be specified if chunk is not None")

        chunk_length = math.floor(n_total_genes / n_chunks)
        chunk_start = chunk * chunk_length
        chunk_end = min(n_total_genes, chunk_start + chunk_length)
        if chunk == n_chunks - 1:
            chunk_end = n_total_genes
    else:
        n_genes = n_total_genes
        chunk_start = 0
        chunk_end = n_genes

    genes = range(chunk_start, chunk_end)
    n_genes = len(genes)
    if n_genes == 0:
        logger.info(
            f"Number of chunks is too large. The pipeline will throw an error beacause there are no genes to test"
        )
    logger.info(f"Processing genes in {genes} from {n_total_genes} in total")
    this_gene_ids = [gene_ids[i] for i in genes]

    with h5py.File(genotype_file, "a") as gf:
        with h5py.File(annotation_file, "a") as af:
            build_data_(
                this_gene_ids, G_full, grouped_annotations, annotation_cols, gf, af
            )

    with h5py.File(x_file, "a") as f:
        f.create_dataset("X", data=X)

    with h5py.File(y_file, "a") as f:
        f.create_dataset("y", data=Y)


if __name__ == "__main__":
    cli()
