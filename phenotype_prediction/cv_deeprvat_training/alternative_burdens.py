import copy
import itertools
import logging
import math
import pickle
import os
import sys
from pathlib import Path
from pprint import pprint
from typing import Dict, List, Optional, Tuple, Union

import click
import dask.dataframe as dd
import numpy as np
import pandas as pd
import psutil
import torch
import torch.nn as nn
import statsmodels.api as sm
import yaml
from joblib import Parallel, delayed
from numcodecs import Blosc
from scipy.sparse import coo_matrix
from statsmodels.tools.tools import add_constant
from torch.utils.data import DataLoader, Dataset, Subset
from tqdm import tqdm
import zarr
import re
from deeprvat.data import DenseGTDataset
from seak import scoretest, lrt
from name_mappings import BTYPES_DICT, PLOF_CONSEQUENCES



logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


PLOF_COLS = PLOF_CONSEQUENCES

BTYPES_DICT_rev = {val: key for key, val in BTYPES_DICT.items()}




@click.group()
def cli():
    pass


def make_dataset_(
    config: Dict,
    debug: bool = False,
    data_key="data",
    btype=None,
    samples: Optional[List[int]] = None,
) -> Dataset:
    data_config = config[data_key]
    if btype is not None:
        rare_config = data_config["dataset_config"]["rare_embedding"]["config"]
        logger.info(f"Adatpting data config to match btype {btype}")
        burden_anno = (
            BTYPES_DICT_rev[btype] if btype in list(BTYPES_DICT_rev.keys()) else btype
        )
        rare_config["annotations"] = [burden_anno]
        if btype == "plof":
            rare_config["annotations"] = ["is_plof"]
            rare_config["thresholds"]["is_plof"] = "is_plof == 1"
            data_config["dataset_config"]["annotations"].append("is_plof")
        else:
            rare_config["thresholds"][
                "Consequence_missense_variant"
            ] = "Consequence_missense_variant == 1"
        logger.info(
            f"rare config: {rare_config}, anntotation file: {data_config['dataset_config']['annotation_file']}"
        )
    ds_pickled = data_config.get("pickled", None)
    if ds_pickled is not None and os.path.isfile(ds_pickled):
        logger.info("Loading pickled dataset")
        with open(ds_pickled, "rb") as f:
            ds = pickle.load(f)
    else:
        variant_file = data_config.get(
            "variant_file", f'{data_config["gt_file"][:-3]}_variants.parquet'
        )
        ds = DenseGTDataset(
            data_config["gt_file"],
            variant_file=variant_file,
            split="",
            skip_y_na=False,
            # None,
            # phenotype_prefix=data_config["gt_file"][:-3],
            **copy.deepcopy(data_config["dataset_config"]),
        )

        restrict_samples = config.get("restrict_samples", None)
        if debug:
            logger.info("Debug flag set; Using only 1000 samples")
            ds = Subset(ds, range(1_000))
        elif samples is not None:
            ds = Subset(ds, samples)
        elif restrict_samples is not None:
            ds = Subset(ds, range(restrict_samples))


    return ds


@cli.command()
@click.option("--debug", is_flag=True)
@click.option("--data-key", type=str, default="data")
@click.option("--btype", type=str)
@click.argument("config-file", type=click.Path(exists=True))
@click.argument("out-file", type=click.Path())
def make_dataset(
    debug: bool,
    data_key: str,  # sample_file: Optional[str],
    config_file: str,
    out_file: str,
    btype: str = None,
):
    with open(config_file) as f:
        config = yaml.safe_load(f)

    # if sample_file is not None:
    #     with open(sample_file, 'rb') as f:
    #         samples = pickle.load(f)['association_samples']
    # else:
    #     samples = None

    # ds = make_dataset_(config, debug=debug, data_key=data_key, samples=samples)

    ds = make_dataset_(config, debug=debug, data_key=data_key, btype=btype)

    with open(out_file, "wb") as f:
        pickle.dump(ds, f)









def get_plof_burden(
    batch: Dict, plof_idx: np.ndarray, btype: str
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    genotypes = np.array(batch["rare_variant_annotations"][:, :, plof_idx, :])
    logger.info(genotypes.shape)
    # assert set(np.unique(genotypes)) == set([0,1])
    logger.info(f"genotypes shape: {genotypes.shape}")
    burdens = np.sum(
        genotypes, axis=2
    )  # get any plof variant by summing across annotation
    logger.info(f"burdens after summing shape: {burdens.shape}")
    if btype == "plof":
        burdens[burdens > 0] = 1
    burdens = np.max(burdens, axis=2)  # sum plof variants for each gene
    logger.info(f"burdens after maxing shape: {burdens.shape}")
    y = batch["y"]
    x = batch["x_phenotypes"]
    return burdens, y, x


# TODO: This is almost the same as compute_burdens_, should be refactored
def compute_plof_burdens_(
    debug: bool,
    config: Dict,
    ds: torch.utils.data.Dataset,
    cache_dir: str,
    btype: str,
    n_chunks: Optional[int] = None,
    chunk: Optional[int] = None,
) -> Tuple[np.ndarray, zarr.core.Array, zarr.core.Array, zarr.core.Array]:
    try:
        data_config = config["plof_data"]
    except:
        data_config = config["data"]

    ds_full = ds.dataset if isinstance(ds, Subset) else ds
    collate_fn = getattr(ds_full, "collate_fn", None)
    n_total_samples = len(ds)
    # import ipdb; ipdb.set_trace()
    if btype == "plof":
        # PLOF_CONSEQUENCES = [f'Consequence_{c}'  for c in ('splice_acceptor_variant', 'splice_donor_variant',
        #             'frameshift_variant', 'stop_gained', 'stop_lost', 'start_lost')]
        logger.info(f"rare annotations: {ds_full.rare_embedding.annotations}")
        if "is_plof" not in ds_full.rare_embedding.annotations:
            plof_idx = [
                ds_full.rare_embedding.annotations.index(c) for c in PLOF_CONSEQUENCES
            ]
        else:
            plof_idx = [ds_full.rare_embedding.annotations.index("is_plof")]
    else:
        annotation_btype = (
            BTYPES_DICT_rev[btype] if btype in list(BTYPES_DICT_rev.keys()) else btype
        )
        plof_idx = [ds_full.rare_embedding.annotations.index(annotation_btype)]
    logger.info(f"plof_idx {plof_idx}")
    if chunk is not None:
        if n_chunks is None:
            raise ValueError("n_chunks must be specified if chunk is not None")

        chunk_length = math.ceil(n_total_samples / n_chunks)
        chunk_start = chunk * chunk_length
        chunk_end = min(n_total_samples, chunk_start + chunk_length)
        samples = range(chunk_start, chunk_end)
        n_samples = len(samples)
        ds = Subset(ds, samples)

        logger.info(f"Processing samples in {samples} from {n_total_samples} in total")
    else:
        n_samples = n_total_samples
        chunk_start = 0
        chunk_end = n_total_samples

    logger.info("Computing burden scores")
    if debug:
        logger.info("Debug mode: Using only 100 samples")
        batch_size = 100
        chunk = chunk if chunk is not None else 0
        n_chunks = n_chunks if n_chunks is not None else 1
        n_total_samples = batch_size * n_chunks
        chunk_start = chunk * batch_size
        chunk_end = chunk_start + batch_size
    else:
        batch_size = n_samples
    print(data_config["dataloader_config"])
    dl = DataLoader(
        ds,
        collate_fn=collate_fn,
        batch_size=batch_size,
        **data_config["dataloader_config"],
    )

    batch = next(iter(dl))
    this_burdens, this_y, this_x = get_plof_burden(batch, plof_idx, btype)

    logger.info("Saving to zarr files")
    burdens = zarr.open(
        Path(cache_dir) / "burdens.zarr",
        mode="a",
        shape=(n_total_samples,) + this_burdens.shape[1:],
        chunks=(1000, 1000),
        dtype=np.float32,
    )

    burdens[chunk_start:chunk_end] = this_burdens
    if this_y.shape[1] > 0:
        # only write results for y and x if y exists.
        # non-existence of y is the case when y_phenotype is empty
        # which can be used in case you just want to get the burdens for all phenotypes
        y = zarr.open(
            Path(cache_dir) / "y.zarr",
            mode="a",
            shape=(n_total_samples,) + this_y.shape[1:],
            chunks=(None, None),
            dtype=np.float32,
        )
        x = zarr.open(
            Path(cache_dir) / "x.zarr",
            mode="a",
            shape=(n_total_samples,) + this_x.shape[1:],
            chunks=(None, None),
            dtype=np.float32,
        )
        # ids = zarr.open(Path(cache_dir) / 'ids.zarr',
        #             mode='a',
        #             shape=(n_total_samples, ) + this_y.shape[1:],
        #             chunks=(None, None),
        #             dtype=np.float32)

        y[chunk_start:chunk_end] = this_y
        x[chunk_start:chunk_end] = this_x
        # ids[chunk_start:chunk_end] = np.expand_dims(np.array([int(i) for i in batch['sample']]),1)

    else:
        y = None
        x = None
    if debug:
        logger.info(
            "Wrote results for chunk indices " f"[{chunk_start}, {chunk_end - 1}]"
        )


    return ds_full.rare_embedding.genes, burdens, y, x


@torch.no_grad()
@cli.command()
@click.option("--debug", is_flag=True)
@click.option("--n-chunks", type=int)
@click.option("--chunk", type=int)
@click.option("--btype", type=str)
@click.option("--dataset-file", type=click.Path(), default=None)
@click.argument("config-file", type=click.Path(exists=True))
@click.argument("out-dir", type=click.Path(exists=True))
def compute_plof_burdens(
    debug: bool,
    n_chunks: Optional[int],
    chunk: Optional[int],
    dataset_file: Optional[str],
    config_file: str,
    out_dir: str,
    btype: str = "plof",
):
    with open(config_file) as f:
        config = yaml.safe_load(f)
    logger.info(f"Extracting burdens for annotation {btype}")
    # ds_rare_config = config['plof_data']['dataset_config']['rare_embedding'][
    #     'config']
    # PLOF_CONSEQUENCES = ds_rare_config['annotations']

    # ds_rare_config['thresholds']['plof'] = ' or '.join(
    #     [f'{c} > 0' for c in PLOF_CONSEQUENCES])
    if dataset_file is not None:
        if Path(dataset_file).is_file():
            logger.info(f"Loading pickled dataset {dataset_file}")
            with open(dataset_file, "rb") as f:
                dataset = pickle.load(f)
    else:
        if dataset_file is None:
            logger.info("No dataset_file provided; instantiating dataset")
        dataset = make_dataset_(config, data_key="plof_data")
        if dataset_file is not None:
            logger.info(f"Pickling dataset to {dataset_file}")
            with open(dataset_file, "wb") as f:
                pickle.dump(dataset, f)

    ds_dataset = dataset.dataset if isinstance(dataset, Subset) else dataset

    genes, _, _, _ = compute_plof_burdens_(
        debug, config, dataset, out_dir, btype, n_chunks=n_chunks, chunk=chunk
    )

    logger.info("Saving computed burdens, corresponding genes, and targets")
    np.save(Path(out_dir) / "genes.npy", genes)










if __name__ == "__main__":
    cli()
