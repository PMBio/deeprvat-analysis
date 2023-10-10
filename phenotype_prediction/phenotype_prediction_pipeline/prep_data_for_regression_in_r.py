import pandas as pd
from pathlib import Path
import zarr
import numpy as np
from scipy.sparse import coo_matrix, spmatrix
import scipy.stats
import statsmodels.api as sm
import pickle
from plotnine import *
from deeprvat.data import DenseGTDataset
import copy
import click
import sys
import os
import logging
import yaml
from prs_comparison import (
    load_prs,
    compare_with_baseline_studies,
    export_all_data,
    get_residual_phenotypes,
    PHENOTYPE_MAPPING,
)
from deeprvat.utils import my_quantile_transform


logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def add_corrected_phenotype_cols(
    phenotype, phenotype_df, prs_df, base_covariate_cols=["age"]
):
    # base_covariate_cols = ['age', *[f'genetic_PC_{i}' for i in range(1,21)]]):
    all_cov_cols = (
        base_covariate_cols
        + ["genetic_sex"]
        + [f"genetic_PC_{i}" for i in range(1, 21)]
    )
    this_phenotype_df = phenotype_df.copy()[
        [*[phenotype, f"{phenotype}_standardized"], *all_cov_cols]
    ].dropna()
    prs_with_pheno = (
        this_phenotype_df.reset_index()
        .rename(columns={"samples": "sample"})
        .set_index("sample")
        .join(prs_df)
        .dropna()
    )

    logger.info("adding quantile transformed phenotype")
    prs_with_pheno[f"{phenotype}_qt"] = my_quantile_transform(prs_with_pheno[phenotype])
    this_phenotype_df[f"{phenotype}_qt"] = my_quantile_transform(
        this_phenotype_df[phenotype]
    )
    res_dict = {}
    logger.info("Computing residual phenotypes")
    for phenotype_suffix in ["", "_standardized", "_qt"]:
        phenotype_col = f"{phenotype}{phenotype_suffix}"
        res_dict[phenotype_col] = {}

        for use_common_prs in [True, False]:
            covariate_cols = copy.deepcopy(base_covariate_cols)

            if phenotype_suffix != "_standardized":
                covariate_cols.append("genetic_sex")

            if use_common_prs:
                covariate_cols.append("common_PRS")

            logger.info(f"Using covariate cols {covariate_cols}")
            Y_hat, residuals, model = get_residual_phenotypes(
                gene_burdens=None,
                ordered_sample_ids=prs_with_pheno.index,
                phenotype_df=prs_with_pheno,
                covariate_cols=covariate_cols,
                phenotype_col=phenotype_col,
                substract_covariates=False,
                fitted_model=None,
            )

            # residuals.hist(alpha = 0.5, bins = 300)
            resid_suffix = "resid_with_prs" if use_common_prs else "resid_wo_prs"
            Y_true = prs_with_pheno[phenotype_col].loc[Y_hat.index]
            residuals.name = f"{phenotype_col}_{resid_suffix}"
            this_phenotype_df = this_phenotype_df.join(residuals)
            dict_key = "with_prs" if use_common_prs else "wo_prs"
            res_dict[phenotype_col][dict_key] = {
                "Y_hat": Y_hat,
                "residuals": residuals,
                "Y_true": Y_true,
            }

    return res_dict, this_phenotype_df


def get_matched_sample_indices(x, y):
    """
    # this function is supposed to do the same as
    # indices= np.array([np.where(x==iy)[0][0] for iy in y]) but is much faster
    #https://stackoverflow.com/questions/8251541/numpy-for-every-element-in-one-array-find-the-index-in-another-array

    Args:
        x : query array
        y: query values. The function returns the index of each element of y in x
    Returns:
        np.array: Index of each element of y in x
    """
    assert np.in1d(y, x).sum() == len(y), "All values of y must be in x"

    xsorted = np.argsort(x)
    ypos = np.searchsorted(x[xsorted], y)
    x_indices = xsorted[ypos]

    x_mask = np.zeros(np.shape(x)).astype(bool)
    x_mask[x_indices] = True

    return x_indices, x_mask


# baseline_indices_to_keep, baseline_mask = get_matched_sample_indices(baseline_samples_int, ordered_samples_to_keep)
def get_sample_indices(dataset_dir, phenotype, split, prs_df):
    logger.info(f"Retrieving index mappings for split {split}")
    split_suffix = "_test" if split == "test" else ""
    with open(f"{dataset_dir}/association_dataset{split_suffix}.pkl", "rb") as f:
        ds_deeprvat = pickle.load(f)
    # with open(f'{dataset_dir_2}/association_dataset{split_suffix}.pkl', 'rb') as f:
    #     ds_deeprvat_2 = pickle.load(f)
    logger.info(
        f"Number of samples before removing na samples: {ds_deeprvat.n_samples}"
    )

    logger.info("Removing samples with missing values if not already done before")
    phenotype_df_cv = ds_deeprvat.phenotype_df
    phenotype_df_cv.index.name = "sample"
    phenotype_df_cv.index = phenotype_df_cv.index.astype(int)
    index_pre_join = phenotype_df_cv.copy().index
    phenotype_df_cv = phenotype_df_cv[phenotype].to_frame().join(prs_df, how="left")
    index_post_join = phenotype_df_cv.copy().index
    assert all(index_pre_join == index_post_join)

    mask_cols = [phenotype, "common_PRS"]
    mask = (phenotype_df_cv[mask_cols].notna()).all(axis=1)
    index_map = np.arange(len(phenotype_df_cv))[mask]

    logger.info(f"Number of samples after removing na samples: {mask.sum()}")

    # check which samples to keep based on sample ids in dense_gt dataset
    # index order corresponds to sample order in deeprvat_burdens and burdens_baseline

    samples_int = [int(i) for i in ds_deeprvat.samples]
    assert all(samples_int == index_pre_join)

    ordered_samples_to_keep = [ds_deeprvat.samples[i] for i in index_map]
    ordered_samples_to_keep = [int(i) for i in ordered_samples_to_keep]

    logger.info(f"Correlation phenotype and PRS: {phenotype_df_cv.corr()}")
    # map to samples in the baseline results. The baseline results comprise all samples and thus
    # need to be subset to only use those in the respective cv split
    # baseline_indices_to_keep, baseline_mask = get_matched_sample_indices(baseline_samples_int,
    #                                                                      ordered_samples_to_keep)
    del ds_deeprvat

    return {"ordered_samples_to_keep": ordered_samples_to_keep, "sample_mask": mask}


@click.group()
def cli():
    pass


def find_substring_position(substring, string_list):
    matches = [index for index, string in enumerate(string_list) if substring in string]
    return matches[0] if matches else -1


@cli.command()
@click.option("--result-dir", type=click.Path(exists=True), default="./")
@click.option("--discoveries-file", type=click.Path(exists=True))
@click.option("--burdens-dir", type=click.Path(exists=True))
@click.option("--dataset-dir", type=click.Path(exists=True))
@click.option("--out-dir", type=click.Path(), default="./postprocessing/r_data")
@click.option("--debug", is_flag=True)
@click.option("--phenotype", type=str)
@click.option("--burden-type", multiple=True)  # , default = ('deeprvat', 'plof'))
@click.argument("config_file", type=click.Path(exists=True))
def prep_r_regression_data(
    result_dir,
    discoveries_file,
    burdens_dir,
    dataset_dir,
    out_dir,
    debug,
    phenotype,
    burden_type,
    config_file,
):
    with open(config_file) as f:
        config = yaml.safe_load(f)

    burden_types = list(burden_type)
    logger.info(f"preparing input data for burden types {burden_types}")
    phenotype = phenotype.replace("_standardized", "")
    logger.info(f"Preparing data for phenotype: {phenotype}")
    prs_mapper = pd.read_csv(config.get("prs_pheno_map"), header=None)

    start_row = find_substring_position("PGS", prs_mapper[1])
    PRS_DICT = dict(zip(prs_mapper.iloc[start_row:, 0], prs_mapper.iloc[start_row:, 1]))

    PRS_DICT = {
        f"{key.replace(' ', '_').replace('-','_')}": value
        for key, value in PRS_DICT.items()
    }
    # PRS_DICT['Total_bilirubin'] = 'PGS001942' #TODO check this
    logger.info(PRS_DICT)
    phenotype_file = config.get("phenotype_file")

    phenotype_df = pd.read_parquet(phenotype_file)
    phenotype_df.index = phenotype_df.index.astype(int)

    splits = ["train", "test"]

    # for phenotype, phenotype_name in phenotypes.items():
    prs_pheno = (
        phenotype
        if phenotype not in list(PHENOTYPE_MAPPING.keys())
        else PHENOTYPE_MAPPING[phenotype]
    )
    print(prs_pheno)
    prs_pheno = prs_pheno.replace(" ", "_").replace("-", "_")
    print(prs_pheno)

    prs_id = PRS_DICT[prs_pheno]
    logger.info(f"Using PRS id {prs_id}")
    prs_df = load_prs(prs_id, config.get("prs_file"))
    logger.info(prs_df)

    ### Define phenotype df
    _, this_phenotype_df = add_corrected_phenotype_cols(phenotype, phenotype_df, prs_df)

    logger.info(this_phenotype_df.columns)

    discoveries = pd.read_parquet(discoveries_file)
    genes_this_pheno = compare_with_baseline_studies(
        phenotype,
        discoveries,
    )

    logger.info("Loading burdens")
    splits = ["train", "test"]
    burdens_dict = {split: {} for split in splits}
    for split in splits:
        for burden_type in burden_types:
            query_burden = burden_type.split("_")[0]
            logger.info(f"Reading burdens for {query_burden}")

            with open(
                f"{burdens_dir}/{split}_{query_burden}_burdens.pickle", "rb"
            ) as f:
                this_burdens = pickle.load(f)
            burdens_dict[split][burden_type] = this_burdens

    indices_dict = {}
    for split in splits:
        indices_dict[split] = get_sample_indices(dataset_dir, phenotype, split, prs_df)
    assert (
        indices_dict["train"]["sample_mask"].shape[0]
        == burdens_dict["train"]["deeprvat"][
            list(burdens_dict["train"]["deeprvat"].keys())[0]
        ][0].shape[0]
    )
    for split in splits:
        data_dict = export_all_data(
            burdens_dict=burdens_dict[split],
            prs_df=prs_df,
            phenotype_df=this_phenotype_df,
            genes_this_pheno=genes_this_pheno,
            ordered_samples_to_keep=indices_dict[split]["ordered_samples_to_keep"],
            sample_mask=indices_dict[split]["sample_mask"],
            covariate_cols=[
                "age",
                "genetic_sex",
                *[f"genetic_PC_{i}" for i in range(1, 21)],
            ],
            phenotype_cols=[
                col for col in this_phenotype_df.columns if phenotype in col
            ],
        )  # the phenotyype stand. for sex group
        # all_data_dict[phenotype_name][fdr][cv_split][split] = data_dict
        logger.info("writing data")
        for method, this_df in data_dict.items():
            for this_key, this_df in this_df.items():
                for gene_list, this_df in this_df.items():
                    for x_y_key, this_df in this_df.items():
                        y_out_file = f"{out_dir}/{split}-y.parquet"
                        this_df.index.name = "sample"
                        print(this_df.head())

                        if x_y_key == "y":
                            # if not Path(y_out_file).exists():
                            #     this_df.to_parquet(y_out_file)
                            this_df.to_parquet(y_out_file)

                        else:
                            # for deeprvat burdens the file path is
                            # {split}-{deeprvat}-{gene_list}-{deeprvat_repeat}-{x_y_key} where gene_list is one of
                            # baseline_discoveries, deeprvat_discoveries (baseline + deeprvat) or genebass
                            # where the output comprises deeprvat burdens for the respective gene list
                            x_out_file = f"{out_dir}/{split}-{method}-{this_key}-{gene_list}-{x_y_key}.parquet"
                            this_df.to_parquet(x_out_file)

        del data_dict

    logger.info("finished")


if __name__ == "__main__":
    cli()