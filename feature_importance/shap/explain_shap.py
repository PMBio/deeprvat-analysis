import copy
import itertools
import logging
import math
import pickle
import os
import sys
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
from scipy.sparse import coo_matrix
from statsmodels.tools.tools import add_constant
from torch.utils.data import DataLoader, Dataset, Subset
from tqdm import tqdm
import zarr
import re
import shap
from matplotlib import pyplot


import deeprvat.deeprvat.models as pl_models


logging.basicConfig(format='[%(asctime)s] %(levelname)s:%(name)s: %(message)s',
                    level='INFO',
                    stream=sys.stdout)
logger = logging.getLogger(__name__)

## SET RANDOM SEED
np.random.seed(0)


class AggregateBagPredictions(pl.LightningDataModule):
    def __init__(self,
                 bag_ckpts):
        
        self.bag_ckpts = bag_ckpts
        
        
    def forward(self, rare_variant_annotations, x_phenotypes):
        preds = sum([ m(rare_variant_annotations, x_phenotypes)
                  for m in self.bag_ckpts]) / len(self.bag_ckpts)
        
        return preds
    
    
#### Multi Pheno Classes    
class MultiphenoDataset(Dataset):
    def __init__(
        self,
        # input_tensor: zarr.core.Array,
        # covariates: zarr.core.Array,
        # y: zarr.core.Array,
        data: Dict[str, Dict],
        min_variant_count: int,
        batch_size: int,
        split: str = "train",
        cache_tensors: bool = False,
        # samples: Optional[Union[slice, np.ndarray]] = None,
        # genes: Optional[Union[slice, np.ndarray]] = None
    ):
        'Initialization'
        super().__init__()

        self.data = data
        self.phenotypes = self.data.keys()
        logger.info(
            f"Initializing MultiphenoDataset with phenotypes:\n{pformat(list(self.phenotypes))}"
        )

        self.cache_tensors = cache_tensors

        for pheno, pheno_data in self.data.items():
            if pheno_data["y"].shape == (
                    pheno_data["input_tensor_zarr"].shape[0], 1):
                pheno_data["y"] = pheno_data["y"].squeeze()
            elif pheno_data["y"].shape != (
                    pheno_data["input_tensor_zarr"].shape[0], ):
                raise NotImplementedError(
                    'Multi-phenotype training is only implemented via multiple y files'
                )

            if self.cache_tensors:
                pheno_data["input_tensor"] = pheno_data["input_tensor_zarr"][:]

        self.min_variant_count = min_variant_count
        self.samples = {
            pheno: pheno_data["samples"][split]
            for pheno, pheno_data in self.data.items()
        }
        self.subset_samples()
        
        print(self.samples)

        self.total_samples = sum([s.shape[0] for s in self.samples.values()])

        self.batch_size = batch_size

        self.sample_order = pd.DataFrame({
            "phenotype":
            itertools.chain(*[[pheno] * len(self.samples[pheno])
                              for pheno in self.phenotypes])
        })
        self.sample_order = self.sample_order.astype(
            {"phenotype": pd.api.types.CategoricalDtype()})
        self.sample_order = self.sample_order.sample(
            n=self.total_samples)  # shuffle
        self.sample_order["index"] = self.sample_order.groupby(
            "phenotype").cumcount()

    def __len__(self):
        'Denotes the total number of batches'
        return math.ceil(len(self.sample_order) / self.batch_size)

    def __getitem__(self, index):
        'Generates one batch of data'

        # TODO:
        # 1. grab min(batch_size, len(self)) (correct this) from computed indices of self.phenotype_order
        # 2. count phenotypes with np.unique
        # 3. return that many samples from that phenotype
        # Figure out how to keep track of how often a phenotype has already occurred!
        # Don't forget to subset by genes!!!!!
        print('inside get item')
        start_idx = index * self.batch_size
        end_idx = min(self.total_samples, start_idx + self.batch_size)
        print(self.total_samples)
        batch_samples = self.sample_order.iloc[start_idx:end_idx]
        samples_by_pheno = batch_samples.groupby("phenotype")

        
        result = dict()
        for pheno, df in samples_by_pheno:
            idx = df["index"].to_numpy()
            print(len(idx))
            annotations = (
                self.data[pheno]["input_tensor"][idx]
                if self.cache_tensors else
                self.data[pheno]["input_tensor_zarr"][idx, :, :, :])

            result[pheno] = {
                "indices": self.samples[pheno][idx],
                "covariates": self.data[pheno]["covariates"][idx],
                "rare_variant_annotations": annotations,
                "y": self.data[pheno]["y"][idx]
            }


        return result

    
    def subset_samples(self):
        for pheno, pheno_data in self.data.items():
            # First sum over annotations (dim2) for each variant in each gene.
            # Then get the number of non-zero values across all variants in all
            # genes for each sample.
            n_samples_orig = self.samples[pheno].shape[0]
            input_tensor = pheno_data["input_tensor_zarr"][
                self.samples[pheno], :, :, :]
            logger.info(input_tensor.shape)
            n_variants_per_sample = np.sum(np.sum(np.array(input_tensor), axis=2) != 0,
                                           axis=(1, 2))
            n_variant_mask = n_variants_per_sample >= self.min_variant_count

            nan_mask = ~pheno_data["y"][self.samples[pheno]].isnan()
            mask = n_variant_mask & nan_mask.numpy()
            self.samples[pheno] = self.samples[pheno][mask]

            logger.info(f"{pheno}: {self.samples[pheno].shape[0]} / "
                        f"{n_samples_orig} samples kept")
    
    
    
class MultiphenoBaggingData(pl.LightningDataModule):
    def __init__(
            self,
            # input_tensor: torch.Tensor,
            # covariates: torch.Tensor,
            # y: torch.Tensor,
            data: Dict[str, Dict],
            train_proportion: float,
            sample_with_replacement: bool = True,
            min_variant_count: int = 1,
            upsampling_factor: int = 1,  # NOTE: Changed
            batch_size: Optional[int] = None,
            train_batch_size: Optional[int] = None, # ADD THIS
            val_batch_size: Optional[int] = None, # ADD THIS
            num_workers: Optional[int] = 0,
            cache_tensors: bool = False):
        logger.info("Intializing datamodule")

        super().__init__()

        if upsampling_factor < 1:
            raise ValueError("upsampling_factor must be at least 1")

        # self.input_tensor = input_tensor
        # self.covariates = covariates
        # self.y = y
        self.train_batch_size = train_batch_size
        self.val_batch_size = val_batch_size
        rng = np.random.default_rng(seed=0)
        
        self.data = data
        #self.n_genes = {
        #    pheno: self.data[pheno]["genes"].shape[0]
        #    for pheno in self.data.keys()
        #}

        # Get the number of annotations and covariates
        # This is the same for all phenotypes, so we can look at the tensors for any one of them
        any_pheno_data = next(iter(self.data.values()))
        self.n_annotations = any_pheno_data["input_tensor_zarr"].shape[2]
        self.n_covariates = any_pheno_data["covariates"].shape[1]

        for pheno, pheno_data in self.data.items():
            n_samples = pheno_data["input_tensor_zarr"].shape[0]
            assert pheno_data["covariates"].shape[0] == n_samples
            assert pheno_data["y"].shape[0] == n_samples
            # self.n_genes = pheno_data["input_tensor_zarr"].shape[1]

            #train_samples = np.arange(0, self.train_batch_size)
            #val_samples = np.arange(self.train_batch_size, n_samples)
            samples = np.arange(n_samples)
            
            n_train_samples = round(n_samples * train_proportion)
            train_samples = np.sort(
                            rng.choice(samples,
                               size=n_train_samples,
                               replace=False))
            
            
            self.data[pheno]["samples"] = {
                    "train": train_samples,
                     "val": np.setdiff1d(samples, train_samples)
                }
                
               
        
        self.save_hyperparameters("min_variant_count", 
                                  "train_batch_size", "val_batch_size", "num_workers", 
                                  "cache_tensors")




    def train_dataloader(self):
        dataset = MultiphenoDataset(
            # self.train_tensor,
            # self.covariates[self.train_samples, :],
            # self.y[self.train_samples, :],
            self.data,
            self.hparams.min_variant_count,
            self.hparams.train_batch_size,
            split="train",
            cache_tensors=self.hparams.cache_tensors)
        return DataLoader(dataset,
                          batch_size=None,
                          num_workers=self.hparams.num_workers)

    def val_dataloader(self):
        dataset = MultiphenoDataset(
            # self.val_tensor,
            # self.covariates[self.val_samples, :],
            # self.y[self.val_samples, :],
            self.data,
            self.hparams.min_variant_count,
            self.hparams.val_batch_size,
            split="val",
            cache_tensors=self.hparams.cache_tensors)
        return DataLoader(dataset,
                          batch_size=None,
                          num_workers=self.hparams.num_workers)
    

@click.group()
def cli():
    pass



@cli.command()
@click.argument('config-file', type=click.Path(exists=True))
@click.argument('checkpoint-files', type=click.Path(exists=True)) #, nargs=-1
@click.argument('input-dir', type=click.Path(exists=True))
@click.argument('repeat-num', type=str)
@click.argument('sampling-batch-no', type=str)
# phenotype_name, input_tensor_file, covariates_file, y_file, training_genes
def feature_importance(
        config_file: str,
        checkpoint_files: str,
        input_dir: str,
        #out_dir: str,
        repeat_num: str,
        sampling_batch_no: str):
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

    complete_models = {}
    logging.info(checkpoint_files)
    logging.info('number of models below:')
    logging.info(len(glob.glob(f'{checkpoint_files}/*.ckpt')))
    dropped = 0
    for ckpt in glob.glob(f'{checkpoint_files}/*.ckpt'):
        if Path(ckpt + '.dropped').is_file():
                # Ignore checkpoints that were chosen to be dropped
            dropped += 1
            continue
        for pheno in config['phenotypes'].keys():
            model = model_class.load_from_checkpoint(
                ckpt,
                config=config['model']['config'],
                pheno=pheno,
                # n_annotations=len(ds_dataset.annotations),
                # n_covariates=len(ds_dataset.x_phenotypes),
                # n_genes=len(ds_dataset.rare_embedding.genes)
            )
            print(model)
            model = model.eval()
            model = model.to(device)
        
            ## TODO: turn this into dictionary 
            complete_models[pheno] = model

        
    logger.info(
            f'Kept {len(complete_models)} models (dropped {dropped}) for repeat {repeat_num}'
        )
    
    ## needs dir to collect input & y from all phenotypes
    compute_shap_values(config,  
              input_dir,
              repeat_num,
              sampling_batch_no,
              complete_models)

    

    
def compute_shap_values(config, 
              input_dir,
              repeat_num,
              sampling_batch_no,
              complete_models):
    out_dir = f'sample_{sampling_batch_no}'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    sampling_batch_no = int(sampling_batch_no)
    samples = slice(None) ## no need for loading from file for now?

    train_batch_size = 3000
    val_batch_size = 1000
    
    shap_covs = {}
    shap_annots = {}
    for pheno in config['phenotypes'].keys():
        
        if os.path.exists(f'{input_dir}/{pheno}'):
            logging.info(f'Loading data for = {pheno}.')
            data = dict()
            data[pheno] = dict()
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

            data[pheno]["training_genes"] = training_genes
            #training_genes_ix = list(training_genes['id'])
            input_tensor = torch.tensor(zarr.open(input_tensor_file,  mode='r')[:])
            #if sampling_batch_no * sample_size < input_tensor.shape[0]:
            data[pheno]["input_tensor_zarr"] =  input_tensor
            data[pheno]["covariates"] = torch.tensor(
                zarr.open(covariates_file, mode='r')[:])
            data[pheno]["y"] = torch.tensor(zarr.open(y_file,
                                                      mode='r')[:])

            logging.info(data[pheno]["input_tensor_zarr"].shape)
            logging.info(data[pheno]["covariates"].shape)
            logging.info(data[pheno]["y"].shape)


            dm = MultiphenoBaggingData( 
                                 data, 
                                 train_proportion = 0.75,
                                 sample_with_replacement = True,
                                 min_variant_count = 1,
                                 upsampling_factor = 1,  
                                 batch_size = None,
                                 train_batch_size = train_batch_size, 
                                 val_batch_size = val_batch_size)

            train_dl = dm.train_dataloader()
            val_dl = dm.val_dataloader()     

            tr_dl_iterator = iter(train_dl)
            for i in range(int(sampling_batch_no)):
                try:
                    background = next(tr_dl_iterator)
                except StopIteration:
                    tr_dl_iterator = iter(train_dl)
                    background = next(tr_dl_iterator)

            val_dl_iterator = iter(val_dl)
            for i in range(int(sampling_batch_no)):
                try:
                    test = next(val_dl_iterator)
                except StopIteration:
                    val_dl_iterator = iter(val_dl)
                    test = next(val_dl_iterator)


            logging.info(background.keys())
            logging.info(f'Training shap explainer, pheno= {pheno}.')
            logging.info(background[pheno]['rare_variant_annotations'].shape)
            e = shap.DeepExplainer(complete_models[pheno], 
                                       [ background[pheno]['rare_variant_annotations'], 
                                         background[pheno]['covariates']] )

            logging.info(f'Generating shap values, pheno= {pheno}.')
            logging.info(test[pheno]['rare_variant_annotations'].shape)
            shap_values = e.shap_values( [test[pheno]['rare_variant_annotations'], 
                                                  test[pheno]['covariates']] )

            logging.info(shap_values[0].shape)
            shap_covs[pheno]=shap_values[1]
            shap_annots[pheno]=shap_values[0]


    with open(f'{out_dir}/{repeat_num}_shap_avg_annots.pkl', 'wb') as f:
        pickle.dump(shap_annots, f)


    with open(f'{out_dir}/{repeat_num}_shap_avg_covs.pkl', 'wb') as ff:
        pickle.dump(shap_covs, ff)


    
    
    
if __name__ == '__main__':
    cli()
   



