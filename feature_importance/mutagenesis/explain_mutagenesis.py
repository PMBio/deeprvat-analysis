import copy
import itertools
import logging
import math
import pickle
import os
import sys
import zarr
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from pprint import pformat, pprint
import glob

import pytorch_lightning as pl
import click
import numpy as np
import pandas as pd
import psutil
import torch
import torch.nn as nn
import statsmodels.api as sm
import yaml

from joblib import Parallel, delayed
from numcodecs import Blosc
from torch.utils.data import DataLoader, Dataset, Subset
from tqdm import tqdm

import deeprvat.deeprvat.models as pl_models


logging.basicConfig(format='[%(asctime)s] %(levelname)s:%(name)s: %(message)s',
                    level='INFO',
                    stream=sys.stdout)
logger = logging.getLogger(__name__)


class MultiphenoDataset(Dataset):
    """
    class used to structure the data and present a __getitem__ function to
    the dataloader, that will be used to load batches into the model
    """

    def __init__(
        self,
        # input_tensor: zarr.core.Array,
        # covariates: zarr.core.Array,
        # y: zarr.core.Array,
        data: Dict[str, Dict],
        # min_variant_count: int,
        batch_size: int,
        split: str = "train",
        cache_tensors: bool = False,
        temp_dir: Optional[str] = None,
        chunksize: int = 1000,
        # samples: Optional[Union[slice, np.ndarray]] = None,
        # genes: Optional[Union[slice, np.ndarray]] = None
    ):
        """
        Initialize the MultiphenoDataset.

        :param data: Underlying dataframe from which data is structured into batches.
        :type data: Dict[str, Dict]
        :param min_variant_count: Minimum number of variants available for each gene.
        :type min_variant_count: int
        :param batch_size: Number of samples/individuals available in one batch.
        :type batch_size: int
        :param split: Contains a prefix indicating the dataset the model operates on. Defaults to "train". (optional)
        :type split: str
        :param cache_tensors: Indicates if samples have been pre-loaded or need to be extracted from zarr. (optional)
        :type cache_tensors: bool
        """

        super().__init__()

        self.data = copy.deepcopy(data)
        self.phenotypes = self.data.keys()
        logger.info(
            f"Initializing MultiphenoDataset with phenotypes:\n{pformat(list(self.phenotypes))}"
        )

        self.cache_tensors = cache_tensors
        if self.cache_tensors:
            self.zarr_root = zarr.group()
        elif temp_dir is not None:
            temp_path = Path(resolve_path_with_env(temp_dir)) / "deeprvat_training"
            temp_path.mkdir(parents=True, exist_ok=True)
            self.input_tensor_dir = TemporaryDirectory(
                prefix="training_data", dir=str(temp_path)
            )
            # Create root group here

        self.chunksize = chunksize
        if self.cache_tensors:
            logger.info("Keeping all input tensors in main memory")

        for pheno, pheno_data in self.data.items():
            if pheno_data["y"].shape == (pheno_data["input_tensor_zarr"].shape[0], 1):
                pheno_data["y"] = pheno_data["y"].squeeze()
            elif pheno_data["y"].shape != (pheno_data["input_tensor_zarr"].shape[0],):
                raise NotImplementedError(
                    "Multi-phenotype training is only implemented via multiple y files"
                )

            if self.cache_tensors:
                zarr.copy(
                    pheno_data["input_tensor_zarr"],
                    self.zarr_root,
                    name=pheno,
                    chunks=(self.chunksize, None, None, None),
                    compressor=Blosc(clevel=1),
                )
                pheno_data["input_tensor_zarr"] = self.zarr_root[pheno]
                # pheno_data["input_tensor"] = pheno_data["input_tensor_zarr"][:]
            elif temp_dir is not None:
                tensor_path = (
                    Path(self.input_tensor_dir.name) / pheno / "input_tensor.zarr"
                )
                zarr.copy(
                    pheno_data["input_tensor_zarr"],
                    zarr.DirectoryStore(tensor_path),
                    chunks=(self.chunksize, None, None, None),
                    compressor=Blosc(clevel=1),
                )
                pheno_data["input_tensor_zarr"] = zarr.open(tensor_path)

        # self.min_variant_count = min_variant_count
        self.samples = {
            pheno: pheno_data["samples"][split]
            for pheno, pheno_data in self.data.items()
        }

        # self.subset_samples()

        self.total_samples = sum([s.shape[0] for s in self.samples.values()])

        self.batch_size = batch_size
        # index all samples and categorize them by phenotype, such that we
        # get a dataframe repreenting a chain of phenotypes
        self.sample_order = pd.DataFrame(
            {
                "phenotype": itertools.chain(
                    *[[pheno] * len(self.samples[pheno]) for pheno in self.phenotypes]
                )
            }
        )
        self.sample_order = self.sample_order.astype(
            {"phenotype": pd.api.types.CategoricalDtype()}
        )
        self.sample_order = self.sample_order.sample(n=self.total_samples)  # shuffle
        # phenotype specific index; e.g. 7. element total, 2. element for phenotype "Urate"
        self.sample_order["index"] = self.sample_order.groupby("phenotype").cumcount()

    def __len__(self):
        "Denotes the total number of batches"
        return math.ceil(len(self.sample_order) / self.batch_size)

    def __getitem__(self, index):
        "Generates one batch of data"

        # 1. grab min(batch_size, len(self)) from computed indices of self.phenotype_order
        # 2. count phenotypes with np.unique
        # 3. return that many samples from that phenotype

        start_idx = index * self.batch_size
        end_idx = min(self.total_samples, start_idx + self.batch_size)
        batch_samples = self.sample_order.iloc[start_idx:end_idx]
        samples_by_pheno = batch_samples.groupby("phenotype", observed=True)

        result = dict()
        for pheno, df in samples_by_pheno:
            # get phenotype specific sub-index
            idx = df["index"].to_numpy()
            assert np.array_equal(idx, np.arange(idx[0], idx[-1] + 1))
            slice_ = slice(idx[0], idx[-1] + 1)

            # annotations = (
            #     self.data[pheno]["input_tensor"][slice_]
            #     if self.cache_tensors
            #     else self.data[pheno]["input_tensor_zarr"][slice_, :, :, :]
            # )
            annotations = self.data[pheno]["input_tensor_zarr"][slice_, :, :, :]

            result[pheno] = {
                "indices": self.samples[pheno][slice_],
                "covariates": self.data[pheno]["covariates"][slice_],
                "rare_variant_annotations": torch.tensor(annotations),
                "y": self.data[pheno]["y"][slice_],
            }

        return result


    def index_input_tensor_zarr(self, pheno: str, indices: np.ndarray):
        # IMPORTANT!!! Never call this function after self.subset_samples()

        x = self.data[pheno]["input_tensor_zarr"]
        first_idx = indices[0]
        last_idx = indices[-1]
        slice_ = slice(first_idx, last_idx + 1)
        arange = np.arange(first_idx, last_idx + 1)
        z = x[slice_]
        slice_indices = np.nonzero(np.isin(arange, indices))
        return z[slice_indices]


class MultiphenoBaggingData(pl.LightningDataModule):
    """
    Preprocess the underlying dataframe, to then load it into a dataset object
    """

    def __init__(
        self,
        data: Dict[str, Dict],
        train_proportion: float,
        sample_with_replacement: bool = True,
        # min_variant_count: int = 1,
        upsampling_factor: int = 1,
        batch_size: Optional[int] = None,
        num_workers: Optional[int] = 0,
        pin_memory: bool = False,
        cache_tensors: bool = False,
        temp_dir: Optional[str] = None,
        chunksize: int = 1000,
    ):
        """
        Initialize the MultiphenoBaggingData.

        :param data: Underlying dataframe from which data structured into batches.
        :type data: Dict[str, Dict]
        :param train_proportion: Percentage by which data is divided into training/validation split.
        :type train_proportion: float
        :param sample_with_replacement: If True, a sample can be selected multiple times in one epoch. Defaults to True. (optional)
        :type sample_with_replacement: bool
        :param min_variant_count: Minimum number of variants available for each gene. Defaults to 1. (optional)
        :type min_variant_count: int
        :param upsampling_factor: Percentual factor by which to upsample data; >= 1. Defaults to 1. (optional)
        :type upsampling_factor: int
        :param batch_size: Number of samples/individuals available in one batch. Defaults to None. (optional)
        :type batch_size: Optional[int]
        :param num_workers: Number of workers simultaneously putting data into RAM. Defaults to 0. (optional)
        :type num_workers: Optional[int]
        :param cache_tensors: Indicates if samples have been pre-loaded or need to be extracted from zarr. Defaults to False. (optional)
        :type cache_tensors: bool
        """
        logger.info("Intializing datamodule")

        super().__init__()

        if upsampling_factor < 1:
            raise ValueError("upsampling_factor must be at least 1")

        self.data = data
        # MODIFIED. COMMENTED OUT. 
        #self.n_genes = {
        #    pheno: self.data[pheno]["genes"].shape[0] for pheno in self.data.keys()
        #}

        # Get the number of annotations and covariates
        # This is the same for all phenotypes, so we can look at the tensors for any one of them
        any_pheno_data = next(iter(self.data.values()))
        self.n_annotations = any_pheno_data["input_tensor_zarr"].shape[2]
        self.n_covariates = any_pheno_data["covariates"].shape[1]
        

        for _, pheno_data in self.data.items():
            n_samples = pheno_data["input_tensor_zarr"].shape[0]
            assert pheno_data["covariates"].shape[0] == n_samples
            assert pheno_data["y"].shape[0] == n_samples

            # TODO: Rewrite this for multiphenotype data
            self.upsampling_factor = upsampling_factor
            if self.upsampling_factor > 1:
                raise NotImplementedError("Upsampling is not yet implemented")

                logger.info(
                    f"Upsampling data with original sample number: {self.y.shape[0]}"
                )
                samples = self.upsample()
                n_samples = self.samples.shape[0]
                logger.info(f"New sample number: {n_samples}")
            else:
                samples = np.arange(n_samples)

            # Sample self.n_samples * train_proportion samples with replacement
            # for training, use all remaining samples for validation
            if train_proportion == 1.0: ### MODIFIED.
                train_samples = samples #self.samples
                #self.val_samples = samples #self.samples
                pheno_data["samples"] = {
                    "train": train_samples
                }
            else:
                n_train_samples = round(n_samples * train_proportion)
                rng = np.random.default_rng()
                # select training samples from the underlying dataframe
                train_samples = np.sort(
                    rng.choice(
                        samples, size=n_train_samples, replace=sample_with_replacement
                    )
                )
                # samples which are not part of train_samples, but in samples
                # are validation samples.
                pheno_data["samples"] = {
                    "train": train_samples,
                    "val": np.setdiff1d(samples, train_samples),
                }

        self.save_hyperparameters(
            # "min_variant_count",
            "train_proportion",
            "batch_size",
            "num_workers",
            "pin_memory",
            "cache_tensors",
            "temp_dir",
            "chunksize",
        )

    def upsample(self) -> np.ndarray:
        """
        does not work at the moment for multi-phenotype training. Needs some minor changes
        to make it work again
        """
        unique_values = self.y.unique()
        if unique_values.size() != torch.Size([2]):
            raise ValueError(
                "Upsampling is only supported for binary y, "
                f"but y has unique values {unique_values}"
            )

        class_indices = [(self.y == v).nonzero(as_tuple=True)[0] for v in unique_values]
        class_sizes = [idx.shape[0] for idx in class_indices]
        minority_class = 0 if class_sizes[0] < class_sizes[1] else 1
        minority_indices = class_indices[minority_class].detach().numpy()
        rng = np.random.default_rng()
        upsampled_indices = rng.choice(
            minority_indices,
            size=(self.upsampling_factor - 1) * class_sizes[minority_class],
        )
        logger.info(f"Minority class: {unique_values[minority_class]}")
        logger.info(f"Minority class size: {class_sizes[minority_class]}")
        logger.info(f"Increasing minority class size by {upsampled_indices.shape[0]}")

        self.samples = upsampled_indices

    def train_dataloader(self):
        """
        trainning samples have been selected, but to structure them and make them load
        as a batch they are packed in a dataset class, which is then wrapped by a
        dataloading object.
        """
        logger.info(
            "Instantiating training dataloader "
            f"with batch size {self.hparams.batch_size}"
        )

        dataset = MultiphenoDataset(
            self.data,
            # self.hparams.min_variant_count,
            self.hparams.batch_size,
            split="train",
            cache_tensors=self.hparams.cache_tensors,
            temp_dir=self.hparams.temp_dir,
            chunksize=self.hparams.chunksize,
        )
        return DataLoader(
            dataset,
            batch_size=None,
            num_workers=self.hparams.num_workers,
            pin_memory=self.hparams.pin_memory,
        )

    def val_dataloader(self):
        """
        validation samples have been selected, but to structure them and make them load
        as a batch they are packed in a dataset class, which is then wrapped by a
        dataloading object.
        """
        logger.info(
            "Instantiating validation dataloader "
            f"with batch size {self.hparams.batch_size}"
        )
        dataset = MultiphenoDataset(
            self.data,
            # self.hparams.min_variant_count,
            self.hparams.batch_size, 
            split="val",
            cache_tensors=self.hparams.cache_tensors,
            temp_dir=self.hparams.temp_dir,
            chunksize=self.hparams.chunksize,
        )
        return DataLoader(
            dataset,
            batch_size=None,
            num_workers=self.hparams.num_workers,
            pin_memory=self.hparams.pin_memory,
        )
    

    
@click.group()
def cli():
    pass



@cli.command()
@click.argument('config-file', type=click.Path(exists=True))
@click.argument('checkpoint-files', type=click.Path(exists=True)) #, nargs=-1
@click.argument('input-dir', type=click.Path(exists=True))
#@click.argument('out-dir', type=click.Path(exists=True))
@click.argument('repeat-num', type=str)
@click.argument('binary-annot-ix', type=str)
def mutagenesis(
        config_file: str,
        checkpoint_files: str,
        input_dir: str,
        repeat_num: str,
        binary_annot_ix: str,
        #training_genes_file:str,
        #training_set
         ):
    if len(checkpoint_files) == 0:
        raise ValueError('At least one checkpoint file must be supplied')

    with open(config_file) as f:
        config = yaml.safe_load(f)

    if torch.cuda.is_available():
        logger.info('Using GPU')
        device = torch.device('cuda')
    else:
        logger.info('Using CPU')
        device = torch.device('cpu')

    logger.info('Loading model and checkpoint')
    model_class = getattr(pl_models, config['model']['type'])

    complete_models, agg_models = [], []
    logging.info(checkpoint_files)
    logging.info(len(glob.glob(f'{checkpoint_files}/*.ckpt')))
    dropped = 0
    for ckpt in glob.glob(f'{checkpoint_files}/*.ckpt'):
        if Path(ckpt + '.dropped').is_file():
                # Ignore checkpoints that were chosen to be dropped
            dropped += 1
            continue

        model = model_class.load_from_checkpoint(
                ckpt,
                config=config['model']['config'],
                # n_annotations=len(ds_dataset.annotations),
                # n_covariates=len(ds_dataset.x_phenotypes),
                # n_genes=len(ds_dataset.rare_embedding.genes)
            )
        model = model.eval()
        model = model.to(device)
        ### TODO extracting only agg model 
        logger.info(model)
        complete_models.append(model)
        #agg_models.append(model.agg_model)

        
    logger.info(
            f'Kept {len(complete_models)} models (dropped {dropped}) for repeat {repeat_num}'
        )

    compute_model_predictions(config,  
                              input_dir,
                              repeat_num,
                              binary_annot_ix,
                              complete_models)

    

    
def compute_model_predictions(config, 
                              input_dir,
                              repeat_num,
                              binary_annot_ix,
                              agg_models):
    out_dir = f'annot_{binary_annot_ix}'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    per_pheno_diff = {p:{} for p in config['training']['phenotypes']}

        
    for pheno in config['training']['phenotypes']:
        
        if os.path.exists(f'{input_dir}/{pheno}'):
            logging.info(f'Loading data for = {pheno}.')
            wildtype_data = dict()
            wildtype_data[pheno] = dict()
            mutated_data = dict()
            mutated_data[pheno] = dict()
            
            input_tensor_file = f'{input_dir}/{pheno}/deeprvat/input_tensor.zarr'
            covariates_file =f'{input_dir}/{pheno}/deeprvat/covariates.zarr'
            y_file = f'{input_dir}/{pheno}/deeprvat/y.zarr'
            training_gene_file = f'{input_dir}/{pheno}/deeprvat/seed_genes.parquet'

            if training_gene_file is not None:
                training_genes = pd.read_parquet(training_gene_file, 
                                                 engine='pyarrow')
                #logging.info(training_genes)
            else:
                training_genes = None

            #training_genes_ix = list(training_genes['id'])
            wildtype_data[pheno]["training_genes"] = training_genes
            mutated_data[pheno]["training_genes"] = training_genes

            input_tensor = torch.tensor(zarr.open(input_tensor_file,  mode='r')[:])
            input_covariates = torch.tensor(zarr.open(covariates_file, mode='r')[:])
            input_y = torch.tensor(zarr.open(y_file, mode='r')[:])
            #if sampling_batch_no * sample_size < input_tensor.shape[0]:

            ######### Work in progress #################
            binary_annot_ix = int(binary_annot_ix)
            per_gene_diff = {f'gene_{i}':[]  
                             for i in range(input_tensor.shape[1])}
            for gene_ix in range(input_tensor.shape[1]):
                ## using samples that has a_i ==1
                sample_ix = (input_tensor[:, gene_ix, binary_annot_ix, :] == 1
                                                                ).sum(dim=1) > 0
                logging.info(f'All #samples {len(sample_ix)}')
                logging.info(f'Gene {gene_ix}: #samples {sum(sample_ix)}')

                if sum(sample_ix) > 1:
                    ## wildtype input for each gene
                    ## each time only one gene changes!
                    wildtype_gene_left = input_tensor[sample_ix, 0:gene_ix, :, :]
                    wildtype_gene_right = input_tensor[sample_ix, gene_ix+1:, :, :]

                    input_wild_type = input_tensor[sample_ix, gene_ix, :, :]
                    input_wild_type = torch.unsqueeze(input_wild_type, dim=1)

                    covarites_wild_type = input_covariates[sample_ix, :]
                    y_wild_type = input_y[sample_ix]

                    #logging.info(f'Wild type: {input_wild_type.shape}')
                    wildtype_data[pheno]["input_tensor_zarr"] = input_tensor[
                                                                sample_ix, :, :, :] 
                    wildtype_data[pheno]["covariates"] = covarites_wild_type
                    wildtype_data[pheno]["y"] = y_wild_type

                    ## mutated input for each gene
                    mutated_left = input_wild_type[:, :, 0:binary_annot_ix, :]
                    mutated_right = input_wild_type[:, :, binary_annot_ix+1:, :]
                    mutated_annot = input_wild_type[:, :, binary_annot_ix, :].unsqueeze(dim=2)
                    mutated_annot[mutated_annot == 1, ] = 0
                    assert mutated_annot.sum() == 0, 'problem with mutagenesis'

                    input_mutated = torch.cat((mutated_left, mutated_annot, mutated_right), dim=2)
                    logging.info(f'Mutated: {input_mutated.shape}')
                    complete_mutated = torch.cat((wildtype_gene_left, input_mutated,
                                                          wildtype_gene_right), dim=1)

                    mutated_data[pheno]["input_tensor_zarr"] =  complete_mutated
                    mutated_data[pheno]["covariates"] = covarites_wild_type
                    mutated_data[pheno]["y"] = y_wild_type


                    dm_wt = MultiphenoBaggingData( 
                             wildtype_data, 
                             train_proportion = 1, 
                             sample_with_replacement = True,
                             #min_variant_count = 1,
                             upsampling_factor = 1,  
                             batch_size = 128
                            )

                    dm_mut = MultiphenoBaggingData( 
                             mutated_data, 
                             train_proportion = 1,
                             sample_with_replacement = True,
                             #min_variant_count = 1,
                             upsampling_factor = 1,  
                             batch_size = 128
                            )


                    ##### do model prediction 
                    val_dl_wt = dm_wt.train_dataloader()   
                    val_dl_mt = dm_mut.train_dataloader()   

                    for i, bag in enumerate(agg_models):
                        logging.info(f'Predictions for bag= {i}.')
                        wt_preds_by_bag, mt_preds_by_bag = [], []
                        #with torch.no_grad():
                        for test in val_dl_wt:
                            predictions = bag(test)[pheno] #dim=2
                            predictions = predictions.detach().numpy() 
                            wt_preds_by_bag.append(predictions)

                        for test in val_dl_mt:
                            predictions = bag(test)[pheno] #dim=2
                            predictions = predictions.detach().numpy() 
                            mt_preds_by_bag.append(predictions)

                        #wt_preds_bag = np.asarray(wt_preds_by_bag, dtype=np.float32)
                        wt_preds_bag = np.concatenate(wt_preds_by_bag)
                        mt_preds_bag = np.concatenate(mt_preds_by_bag) #vstack
                        logging.info(f'complete mutant_preds: {mt_preds_bag.shape}')
                        logging.info(f'complete wildty_preds: {wt_preds_bag.shape}')
                        per_sample_diff = np.absolute(wt_preds_bag - mt_preds_bag).squeeze()

                        mean_sample_diff = sum(per_sample_diff) / len(per_sample_diff)
                        per_gene_diff[f'gene_{gene_ix}'].append(mean_sample_diff)

            per_pheno_diff[pheno] = per_gene_diff

    pheno_results = dict()
    for pheno in per_pheno_diff.keys():
        per_gene = per_pheno_diff[pheno]
        gene_avg = []
        for g in per_gene.keys():
            if len(per_gene[g]) >= 1:
                gene_avg.append(per_gene[g])
        
        gene_avg = np.array(gene_avg)
        if len(gene_avg) > 0:
            pheno_results[pheno] = sum(gene_avg) / len(gene_avg)
        else:
             pheno_results[pheno] = [0]
    
    logging.info(per_pheno_diff)
    logging.info(pheno_results)
    with open(f'{out_dir}/diff_{repeat_num}.pkl', 'wb') as f:
        pickle.dump(pheno_results, f)


    
if __name__ == '__main__':
    cli()
   
    
   
