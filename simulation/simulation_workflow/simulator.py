import copy
import pickle
import logging
import os
import sys
from typing import Any, Dict, List, Optional, Set, Tuple, Union
import yaml
import random
from pathlib import Path
from scipy.stats import halfnorm, norm, bernoulli, dirichlet
from scipy.special import logit
import torch.nn.functional as F


import dask.dataframe as dd

import numpy as np
import pandas as pd
import torch
import zarr
from numcodecs import Blosc
from torch.utils.data import DataLoader, Subset
from tqdm import tqdm

import sim_agg_models as sim_agg_models

import deeprvat.data.rare as rare_embedders
from deeprvat.utils import standardize_series
from deeprvat.data import DenseGTDataset

from torch.utils.data import DataLoader


logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

# DEFAULT_CHROMOSOMES = [f'chr{x}' for x in range(1, 23)]

GENE_WEIGHTS = {"normal": norm, "half_normal": halfnorm, "dirichlet": dirichlet}

NOISE = {"normal": norm}

COVARIATE = {"normal": norm, "bernoulli": bernoulli}

VARIANT_WEIGHTS = {"normal": norm, "half_normal": halfnorm, "dirichlet": dirichlet}


class SimulatedPhenotypes:
    def __init__(
        self,
        causal_gene_file: str = None,
        input_tensor_file: str = None,
        standardize_sim_phenotypes: Optional[bool] = True,
        annos_for_causal_var_selection: Optional[List[str]] = None,
        annotation_weight: str = "inverse_normalized",
        prop_causal_variants: int = 0.15,
        max_rare_af: float = 0.01,
        var_expl_variant_noise: float = 0.1,
        var_expl_maf: float = 0.3,
        var_expl_binary: float = 0.3,
        var_expl_cont: float = 0.3,
        var_extra_group: float = 0.0,
        binary_weights: dict = None,
        binary_anno_weight_exp: int = 1,
        # variant_effect_size: float = 0.5,
        var_genetic_effects: float = 0.4,
        var_noise: float = 0.6,
        n_samples: Optional[int] = None,
        gene_weight: Optional[Dict] = None,
        n_causal_genes: Optional[int] = None,
        noise: Optional[dict] = None,
        covariate: Optional[dict] = None,
        sample_variants_only: bool = False,
        config: dict = None,
        debug: bool = False,
    ):
        # select how many samples to use in the simulation
        # something like use annotations (and select which ones)
        # something to define how to compute the phenotypes
        # files/base data frames
        # for large data frames, we provide the option to provide the data frames directly
        # (so they don't need to be loaded when generating this object but have been loaded
        # previously)
        self.debug = debug
        self.sample_variants_only = sample_variants_only
        self.n_samples = n_samples

        # configs
        try:
            self.data_config = config["data"]
        except:
            self.data_config = config["simulation_data"]
        self.simulation_config = config["simulation_config"]

        self.input_tensor_file = input_tensor_file
        # self.input_tensor_file = None

        # gene related
        self.n_causal_genes = n_causal_genes
        if causal_gene_file is not None:
            logger.info("loading causal genes")
            self.prespec_causal_genes = pd.read_parquet(causal_gene_file)
        else:
            raise ValueError(f"Causal_gene file must be defined")

        # variant related
        self.af_col = list(self.data_config["dataset_config"]["min_common_af"].keys())[
            0
        ]

        self.max_rare_af = max_rare_af
        self.var_expl_variant_noise = var_expl_variant_noise
        self.var_expl_binary = var_expl_binary
        self.binary_weights = binary_weights
        self.binary_anno_weight_exp = binary_anno_weight_exp
        self.var_expl_maf = var_expl_maf
        self.var_expl_cont = var_expl_cont
        self.var_extra_group = 0 if var_extra_group is None else var_extra_group

        self.extra_group_annotations = self.simulation_config.get(
            "extra_group_annotations", None
        )
        if (self.var_extra_group > 0) and (self.extra_group_annotations is None):
            ValueError(
                'If var_extra_group is non-zero "extra_group_annotations" must be defined!'
            )

        self.prop_causal_variants = prop_causal_variants
        # self.variant_effect_size = variant_effect_size

        ## variance explained by different components

        self.var_genetic_effects = var_genetic_effects
        self.var_non_genetic_effecs = 1 - var_genetic_effects
        self.var_noise = self.var_non_genetic_effecs * var_noise
        self.var_covariates = self.var_non_genetic_effecs * (1 - var_noise)

        self.annos_for_causal_var_selection = (
            annos_for_causal_var_selection
            if annos_for_causal_var_selection is not None
            else copy.deepcopy(
                self.data_config["dataset_config"]["rare_embedding"]["config"][
                    "annotations"
                ]
            )
        )
        if (self.extra_group_annotations is not None) and (self.var_extra_group > 0):
            self.annos_for_causal_var_selection = (
                self.annos_for_causal_var_selection + self.extra_group_annotations
            )
            # self.all_annotations = set(self.data_config['dataset_config']['annotations']).union(set(self.annos_for_causal_var_selection))
            logger.info(f"Updating config to inlcude additonal annotations. ")
            self.data_config["dataset_config"][
                "annotations"
            ] = self.annos_for_causal_var_selection.copy()
            self.data_config["dataset_config"]["rare_embedding"]["config"][
                "annotations"
            ] = self.annos_for_causal_var_selection.copy()
            print(self.data_config["dataset_config"]["annotations"])

        self.annotation_weight_config = annotation_weight

        logger.info(
            f"Anntoations used for causal variant selection: \
             {self.annos_for_causal_var_selection}"
        )

        self.standardize_sim_phenotypes = standardize_sim_phenotypes

        # noise and distributions:
        self.gene_weight_config = convert_argument(gene_weight)
        self.noise_config = convert_argument(noise)
        self.covariate_config = convert_argument(covariate)

        self.setup_vars_to_ensgid_mapping(
            self.simulation_config["vars_to_ensgid_mapping_file"]
        )

        # genes_to_keep = genes which are in the variant data frame
        logger.info("setting up annotations")
        self.setup_annotations(self.data_config["dataset_config"]["annotation_file"])

    #### Methods for setting up the data #################################################
    ########################################################################################

    def setup_vars_to_ensgid_mapping(self, vars_to_ensgid_mapping_filename):
        logger.info("Reading variant to gene id mapping data frame")
        self.vars_to_ensgid_mapping = pd.read_parquet(
            vars_to_ensgid_mapping_filename, engine="pyarrow"
        )
        self.vars_to_ensgid_mapping = self.vars_to_ensgid_mapping.astype(
            {"variant_id": "int64"}
        ).drop(columns="gene_base")

    def setup_annotations(self, annotation_filename):
        columns = set(["id"] + [self.af_col] + self.annos_for_causal_var_selection)

        logger.info("Reading annotation data frame")
        self.annotation_df = dd.read_parquet(
            annotation_filename, columns=list(set(columns)), engine="pyarrow"
        ).compute()
        self.annotation_df = self.annotation_df.set_index("id")

        annotations_to_negate = ["sift_score"]
        for col in annotations_to_negate:
            if col in self.annotation_df.columns:
                self.annotation_df[col] = -self.annotation_df[col]

        self.binary_annotations = [
            col
            for col in self.annotation_df
            if np.isin(self.annotation_df[col].unique(), [0, 1]).all()
        ]
        logger.info(f"Binary annotation columns: {self.binary_annotations}")
        anno_frequency = self.annotation_df[self.binary_annotations].sum(axis=0)

        logger.info(
            f"Exponent used to weight binary annotations: anno_freq^(-1/{self.binary_anno_weight_exp})"
        )
        self.inverse_anno_weights = anno_frequency ** (-1 / self.binary_anno_weight_exp)
        self.normalized_inverse_anno_weights = self.inverse_anno_weights / (
            (self.inverse_anno_weights).sum()
        )
        # normalized here means that all weights sum up to 1

    #### Methods for simulating Phenotypes #################################################
    ########################################################################################

    # wrapper function
    def simulate_data_set(self):
        logger.info("Simulating data set")

        logger.info("Selecting causal genes and variants")

        # TODO: implement random causal gene selection
        logger.info("Using prespecified causal genes")
        self.causal_genes = self.prespec_causal_genes

        self.vars_causal_genes = self.get_variants_for_causal_genes(
            self.causal_genes["gene"]
        )
        # n_causal_vars_per_gene = self.get_number_of_causal_vars_per_gene(self.vars_causal_genes)

        logger.info("Sampling causal variants")
        self.all_causal_var_ids, self.n_causal_vars_per_gene = self.sample_causal_vars(
            self.vars_causal_genes
        )

        if self.sample_variants_only:
            # option to only analyze the simulated causal variants
            logger.info("Stopping after sampling variants")

        else:
            logger.info(
                f"Total number of causal variants: {len(self.all_causal_var_ids)}"
            )
            logger.info("Embedding variants")

            (
                self.rare_variant_embedding,
                self.samples,
            ) = self.setup_rare_variant_embedding(self.all_causal_var_ids)

            logger.info("Computing gene scores")
            self.gene_burden_scores = self.compute_gene_score(
                self.rare_variant_embedding
            )

            self.sim_phenotypes = self.compute_phenotype(self.gene_burden_scores)

            # Recording of gene weights is buggy
            # self.causal_genes['weight'] = self.gene_weights
            self.phenotype_df = pd.DataFrame(
                {
                    "sim_phenotype": self.sim_phenotypes["sim_phenotype"],
                    **self.covariates,
                },
                index=self.samples["sample_id"],
            )

    def get_variants_for_causal_genes(self, causal_gene_names):
        self.vars_causal_genes = self.vars_to_ensgid_mapping[
            self.vars_to_ensgid_mapping["gene"].isin(causal_gene_names)
        ]
        logger.info(
            f"Number of variants before removing duplicated variants: {len(self.vars_causal_genes)}"
        )

        #########################################
        # remove variants that are in >1 gene (not restricted to simulated causal genes but any gene)
        unique_vars = (
            self.vars_to_ensgid_mapping[
                self.vars_to_ensgid_mapping["variant_id"].isin(
                    self.vars_causal_genes["variant_id"]
                )
            ]
            .groupby("variant_id")
            .size()
            .reset_index(name="n_genes")
        )

        unique_vars = unique_vars[unique_vars["n_genes"] == 1]["variant_id"]

        self.vars_causal_genes = self.vars_causal_genes[
            self.vars_causal_genes["variant_id"].isin(unique_vars)
        ]
        logger.info(
            f"Number of variants after removing duplicated variants: {len(self.vars_causal_genes)}"
        )
        #########################################

        self.annotation_df = self.annotation_df.loc[
            self.vars_causal_genes["variant_id"]
        ]
        # Filter variants by MAF (if you don't want to do any MAF filtering put cutoff to 1)
        # this is the mask for identifying common varians based on MAF as used in DenseGTDataset
        # rare variants are the inverse of the mask
        # af_mask = ((self.annotation_df[self.af_col] >= self.max_rare_af) &
        #             (self.annotation_df[self.af_col] <= 1 - self.max_rare_af))
        af_mask = (self.annotation_df[self.af_col] < self.max_rare_af) & (
            self.annotation_df[self.af_col] > 0
        )

        rare_variants = set(self.annotation_df[af_mask].index)
        n_rare_variants = len(rare_variants)

        logger.info(
            f"Number of variants after AF filtering with MAF < {self.max_rare_af}: {n_rare_variants}"
        )
        self.n_total_causal_variants = round(
            n_rare_variants * self.prop_causal_variants
        )
        logger.info(
            f"Using {self.prop_causal_variants} from all {n_rare_variants} as causal variants, "
            f"i.e., {self.n_total_causal_variants} causal variants"
        )

        self.vars_causal_genes = self.vars_causal_genes[
            self.vars_causal_genes["variant_id"].isin(rare_variants)
        ]

        return self.vars_causal_genes

    def sample_causal_vars(self, vars_causal_genes):
        def rescale_variance_np(X, target_variance):
            # This function is inspired by the rescaleVariance function
            # from the PhenotypeSimulator package (Meyer et al., 2018).
            # It scales the variance such that the mean variance of all columns of
            # X is equal to the target_variance.
            # X can be a one or two dimensional array (but we only have 1d input here)
            mean_variance = np.var(X, axis=0).mean()
            scale_factor = np.sqrt(target_variance / mean_variance)
            X_scaled = scale_factor * X
            return X_scaled, scale_factor

        #################################
        # Transform MAF related column
        # TODO: also implement other transformation such as Beta(1,25)
        maf_col = [col for col in self.annos_for_causal_var_selection if "AF" in col][0]
        logger.info(f"Taking log10 of maf_col: {maf_col}")
        ##
        self.annotation_df[f"{maf_col}_log10"] = -np.log10(self.annotation_df[maf_col])

        self.annos_for_causal_var_selection.remove(maf_col)

        cont_annotations = [
            col
            for col in self.annos_for_causal_var_selection
            if col not in set(self.binary_annotations).union(set([maf_col]))
        ]

        self.annos_for_causal_var_selection = [
            f"{maf_col}_log10"
        ] + self.annos_for_causal_var_selection
        logger.info(
            f"annos_for_causal_var_selection: {self.annos_for_causal_var_selection}"
        )

        #################################
        # score variants

        # var_expl_not_noise = 1 - self.var_expl_variant_noise
        # self.var_expl_maf = round(self.var_expl_maf * var_expl_not_noise, 3)
        # self.var_expl_binary = round((1 - self.var_expl_maf -  self.var_expl_variant_noise) * self.var_expl_binary, 3)
        # self.var_expl_cont = round(1 - self.var_expl_maf -  self.var_expl_variant_noise - self.var_expl_binary, 3)

        variant_variances = {
            "var_expl_variant_noise": self.var_expl_variant_noise,
            "var_expl_maf": self.var_expl_maf,
            "var_expl_binary": self.var_expl_binary,
            "var_expl_cont": self.var_expl_cont,
            "var_extra_group": self.var_extra_group,
        }
        logger.info(
            f"variant_variances: {variant_variances} \n, sum of variances: {np.array(list(variant_variances.values())).sum()}"
        )
        assert round(np.array(list(variant_variances.values())).sum(), 1) == 1

        # split variance across all annotations in the different categories
        if self.binary_weights is not None:
            bin_anno_variance = self.binary_weights
        else:
            bin_anno_variance = (
                self.normalized_inverse_anno_weights * self.var_expl_binary
            ).to_dict()

        if self.extra_group_annotations is not None:
            cont_annotations = [
                col
                for col in cont_annotations
                if col not in self.extra_group_annotations
            ]
            logger.info(
                f"cont annotations after removing extra_group_annotations: {cont_annotations}"
            )

        if self.var_extra_group > 0:
            extra_group_annotation_variance = (
                np.ones(len(self.extra_group_annotations))
                / len(self.extra_group_annotations)
                * self.var_extra_group
            )
            extra_group_annotation_variance = dict(
                zip(self.extra_group_annotations, extra_group_annotation_variance)
            )

        cont_anno_variance = (
            np.ones(len(cont_annotations)) / len(cont_annotations) * self.var_expl_cont
        )
        cont_anno_variance = dict(zip(cont_annotations, cont_anno_variance))

        self.anno_weights = {"UKB_AF_log10": self.var_expl_maf}
        self.anno_weights.update(bin_anno_variance)
        self.anno_weights.update(cont_anno_variance)
        if self.var_extra_group > 0:
            self.anno_weights.update(extra_group_annotation_variance)

        vars_with_anno = self.annotation_df[self.anno_weights.keys()].loc[
            self.vars_causal_genes["variant_id"]
        ]
        logger.info(vars_with_anno.head())
        # rescale variances

        self.vars_with_anno_scaled = vars_with_anno.copy()
        self.scale_factors = {}
        for col in vars_with_anno.columns:
            self.vars_with_anno_scaled[col], scale_factor = rescale_variance_np(
                self.vars_with_anno_scaled[col], self.anno_weights[col]
            )
            self.scale_factors[col] = scale_factor

        self.vars_with_anno_scaled["noise"] = np.random.normal(
            0, self.var_expl_variant_noise, len(self.vars_with_anno_scaled)
        )
        self.var_scores = self.vars_with_anno_scaled.sum(axis=1)
        logger.info(
            f"Number of variants with var scores < 0: {(self.var_scores < 0).sum()}"
        )

        logger.info(f"Variance variant scores: {np.var(self.var_scores)}")

        # select top ranked variants as causal variants
        all_causal_var_ids = (
            self.var_scores.sort_values(ascending=False)
            .head(self.n_total_causal_variants)
            .index.to_list()
        )

        # assign selected causal variants to genes
        causal_vars_per_gene = {}
        genes_wo_vars = []
        for gene in self.vars_causal_genes["gene"].unique():
            print(gene)
            causal_vars_this_gene = self.vars_causal_genes.query(
                "gene == @gene and variant_id in @all_causal_var_ids"
            )["variant_id"]
            if len(causal_vars_this_gene) == 0:
                genes_wo_vars.append(gene)
            else:
                causal_vars_per_gene[gene] = list(causal_vars_this_gene)
        n_causal_vars_per_gene = {
            gene: len(causal_vars_per_gene[gene])
            for gene in causal_vars_per_gene.keys()
        }
        logger.info(f"Number of simulated causal variants: {len(all_causal_var_ids)}")
        logger.info(f"Number of causal variants per gene: {n_causal_vars_per_gene}")
        logger.info(
            f"Genes not carrying any simulated causal variant: {genes_wo_vars}, Removing these genes from causal genes"
        )
        self.causal_genes = self.causal_genes[
            ~self.causal_genes["gene"].isin(genes_wo_vars)
        ]

        return all_causal_var_ids, n_causal_vars_per_gene

    def setup_rare_variant_embedding(self, all_causal_var_ids):
        if self.input_tensor_file is not None and Path(self.input_tensor_file).exists():
            logger.info("Loading saved input tensor")
            input_tensor = torch.tensor(zarr.open(self.input_tensor_file, mode="r")[:])
            samples_path = f"{os.path.dirname(self.input_tensor_file)}/samples.parquet"
            samples = pd.read_parquet(samples_path)
            # try:
            #     samples_path = f'{os.path.dirname(self.input_tensor_file)}/samples.parquet'
            #     samples = pd.read_parquet(samples_path)
            # except: #TODO: remove this in futre (did renaming but actually the pickle file is a parquet file)
            #     samples_path = f'{os.path.dirname(self.input_tensor_file)}/samples.pickle'
            #     samples = pd.read_parquet(samples_path)
        else:
            logger.info("Generating DenseGT Data set")
            self.data_config["dataset_config"]["x_phenotypes"] = ["genetic_sex"]
            self.ds = DenseGTDataset(
                gt_file=self.data_config["gt_file"],
                variant_file=self.data_config["variant_file"],
                variants_to_keep=all_causal_var_ids,
                split="",
                skip_y_na=False,
                skip_x_na=False,
                **self.data_config["dataset_config"],
            )

            collate_fn = self.ds.collate_fn
            pad_value = self.ds.rare_embedding.pad_value

            if self.debug:
                logger.info("Debug flag set; Using only 1000 samples")
                self.ds = Subset(self.ds, range(1_000))
            elif self.n_samples is not None:
                logger.info(f"n_samples set; Using only {self.n_samples} samples")
                self.ds = Subset(self.ds, range(self.n_samples))

            dl = DataLoader(
                self.ds, collate_fn=collate_fn, **self.data_config["dataloader_config"]
            )

            logger.info("Generating input tensor")
            batches = [
                batch
                for batch in tqdm(
                    dl,
                    file=sys.stdout,
                    total=len(self.ds)
                    // self.data_config["dataloader_config"]["batch_size"],
                )
            ]
            self.rare_batches = [b["rare_variant_annotations"] for b in batches]
            samples = [b["sample"] for b in batches]
            samples = [
                index for sublist in [b["sample"] for b in batches] for index in sublist
            ]
            # we only retrieve genetic sex as covariates
            sex = [
                index.item()
                for sublist in [b["x_phenotypes"] for b in batches]
                for index in sublist
            ]
            sex = pd.Series(sex, index=samples, name="genetic_sex").to_frame()
            samples = sex.reset_index().rename(columns={"index": "sample_id"})

            max_n_variants = max(r.shape[-1] for r in self.rare_batches)
            logger.info("Building input_tensor")
            input_tensor = torch.cat(
                [
                    F.pad(r, (0, max_n_variants - r.shape[-1]), value=pad_value)
                    for r in tqdm(self.rare_batches, file=sys.stdout)
                ]
            )

            if self.input_tensor_file is not None:
                logger.info(f"saving input tensor to {self.input_tensor_file}")
                zarr.save_array(
                    self.input_tensor_file,
                    input_tensor.numpy(),
                    chunks=(1000, None, None, None),
                    compressor=Blosc(clevel=1),
                )
                samples_path = (
                    f"{os.path.dirname(self.input_tensor_file)}/samples.parquet"
                )
                samples.to_parquet(samples_path)

        logger.info(f"Input tensor shape {input_tensor.shape} ")

        return input_tensor, samples

    def compute_gene_score(self, rare_variant_embedding):
        model_config = self.simulation_config["sim_agg_model"]
        logger.info("check if these are in the same order!")
        logger.info(self.annos_for_causal_var_selection)
        logger.info(
            self.data_config["dataset_config"]["rare_embedding"]["config"][
                "annotations"
            ]
        )
        weights = np.array(
            [self.scale_factors[col] for col in self.annos_for_causal_var_selection]
        )
        logger.info("computing gene score")
        logger.info(weights)
        logger.info(model_config)
        aggregation_model = getattr(sim_agg_models, model_config["type"])(
            weights=weights, **model_config["config"]
        )

        # gene_burden_scores = aggregation_model(rare_variant_embedding)

        all_gene_burdens = []
        rare_batches = torch.split(rare_variant_embedding, 64, 0)
        logger.info(
            f"Number of rare batches: {len(rare_batches)}, batch zero shape: {rare_batches[0].shape}"
        )
        for batch in tqdm(rare_batches, file=sys.stdout, total=len(rare_batches)):
            this_burden = aggregation_model(batch)
            all_gene_burdens.append(this_burden)
        gene_burden_scores = torch.cat(all_gene_burdens)
        logger.info(f"gene_burden_scores shape: {gene_burden_scores.shape}")
        return gene_burden_scores

    def compute_phenotype(self, gene_burden_scores):
        if self.gene_weight_config is not None:
            logger.info(f"Sampling gene-specific weights: {self.gene_weight_config}")
            n_genes = gene_burden_scores.shape[1]
            logger.info(f"Number of genes: {n_genes}")
            if self.gene_weight_config["type"] == "dirichlet":
                logger.info("Sampling from dirichlet with constant weights")
                alpha_multipl = self.gene_weight_config["config"]["alpha_mutipl"]
                alpha = np.full(n_genes, 1 / n_genes) * alpha_multipl
                dirichlet_config = {"type": "dirichlet", "config": {"alpha": alpha}}
                self.gene_weights = sample_from_distribution(
                    dirichlet_config, GENE_WEIGHTS, 1
                )[0]
            else:
                self.gene_weights = sample_from_distribution(
                    self.gene_weight_config, GENE_WEIGHTS, gene_burden_scores.shape[1]
                )
        else:
            logger.info("Not using gene-specific weights")
            self.gene_weights = torch.full(
                torch.Size([gene_burden_scores.shape[1]]), 1
            )  # make a tensor of length 33

        logger.info(f"Number of genes in input tensor: {len(self.gene_weights)}")
        # phenotype_y = w_variant [\sum_k geneweight_k * burden_k] + covariates + noise
        # pheno_per_sample = torch.sum(torch.mul(gene_burden_scores, self.gene_weights), dim = 1)
        genetic_effects = torch.sum(
            torch.mul(gene_burden_scores, self.gene_weights), dim=1
        )
        n_samples = genetic_effects.shape[0]

        ##Scale variance of genetic effects and add scaled covariance and noise
        logger.info(
            f"Resecaling phenotype components. Variance explained genetic effects: \
            {self.var_genetic_effects}, noise: {self.var_noise}, covariates: {self.var_covariates}"
        )

        assert (
            round(
                round(self.var_genetic_effects, 2)
                + round(self.var_noise, 2)
                + round(self.var_covariates, 2)
            )
            == 1
        )
        self.scaled_genetic_effects = rescale_variance(
            genetic_effects, self.var_genetic_effects
        )

        self.scaled_noise = self.get_noise(n_samples)
        self.scaled_convariates = self.get_covariates(n_samples)
        pheno_per_sample = (
            self.scaled_genetic_effects + self.scaled_convariates + self.scaled_noise
        )
        pheno_per_sample = pd.Series(
            pheno_per_sample,
            index=self.samples["sample_id"],
            name="sim_phenotype_non_std",
        )

        if self.standardize_sim_phenotypes:
            logger.info("Standardizing phenotypes by genetic sex group")
            # mean centering and scaling by standard deviation
            # pheno_per_sample = standardize_series(pheno_per_sample)
            # self.samples.set_index('sample_id') is the same as self.sex
            pheno_per_sample = pheno_per_sample.to_frame().join(
                self.samples.set_index("sample_id")
            )
            pheno_per_sample["sim_phenotype"] = pheno_per_sample.groupby(
                ["genetic_sex"]
            )["sim_phenotype_non_std"].apply(lambda x: standardize_series(x))

        logger.info(f"Simulated phenotype values:{pheno_per_sample}")

        return pheno_per_sample

    def get_noise(self, n_samples):
        if self.noise_config is not None:
            logger.info(f"Adding noise: {self.noise_config}")
            noise_values = sample_from_distribution(self.noise_config, NOISE, n_samples)
            noise_values = rescale_variance(noise_values, self.var_noise)

        return noise_values

    def get_covariates(self, n_samples):
        self.covariates = {}
        if self.covariate_config is not None:
            logger.info(f"Adding covariates {self.covariate_config}")

            cov_list = []
            for cov in self.covariate_config.keys():
                cov_values_with_prefac, cov_values = sample_from_distribution(
                    self.covariate_config[cov],
                    COVARIATE,
                    n_samples,
                    return_raw_values=True,
                )
                # TODO: Do I need to rescale the cov_values as well?
                # not really I think: cov_values_with_prefac is basically
                # cov_value_with_prefac = cov_value * var_scaling_facotr * pre_fac_formula
                # where pre_fac_formula is eg. 0.5 in Y = 0.5 X_1 + 0.5_X_2 + genetic_comp + noise
                self.covariates[f"sim_{cov}"] = pd.Series(
                    cov_values, index=self.samples["sample_id"]
                )
                cov_list.append(
                    cov_values_with_prefac.unsqueeze(dim=1)
                )  # n_samples x1 tensor
        cov_values_with_prefac = rescale_variance(
            torch.cat(cov_list, dim=1).sum(dim=1), self.var_covariates
        )

        return cov_values_with_prefac


def sample_from_distribution(config, class_dict, n, return_raw_values=False):
    dist_class = class_dict[config["type"]]
    pre_fac = config.get("pre_fac", 1)
    class_config = config.get("config", {})
    dist_class = dist_class(**class_config)
    dist_values = torch.Tensor(dist_class.rvs(n))
    dist_values_with_prefac = pre_fac * dist_values

    if return_raw_values:
        return dist_values_with_prefac, dist_values
    return dist_values_with_prefac


def rescale_variance(X, target_variance):
    # This function is inspired by the rescaleVariance function
    # from the PhenotypeSimulator package (Meyer et al., 2018).
    # It scales the variance such that the mean variance of all columns of
    # X is equal to the target_variance.
    # X can be a one or two dimensional array (but we only have 1d input here)
    mean_variance = torch.var(X, axis=0).mean()
    scale_factor = torch.sqrt(target_variance / mean_variance)
    X_scaled = scale_factor * X
    return X_scaled


def convert_argument(argument):
    if argument == "None":
        conv_arg = None
    elif type(argument) == str:
        conv_arg = eval(argument)
    else:
        conv_arg = argument

    return conv_arg
