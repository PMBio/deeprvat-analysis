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
    phenotype, phenotype_df, prs_df, config,# base_covariate_cols=["age"]
):
    # base_covariate_cols = ['age', *[f'genetic_PC_{i}' for i in range(1,21)]]):
    all_cov_cols = config['covariates']
    base_covariate_cols = copy.deepcopy(all_cov_cols)

    this_phenotype_df = phenotype_df.copy()[[phenotype, *all_cov_cols]].dropna()
    prs_with_pheno = (
        this_phenotype_df\
        .join(prs_df)
        .dropna()
    )
    
    # comment this back in !!! 
    # logger.info("adding quantile transformed phenotype")
    # prs_with_pheno[f"{phenotype}_qt"] = my_quantile_transform(prs_with_pheno[phenotype])
    # this_phenotype_df[f"{phenotype}_qt"] = my_quantile_transform(
    #     this_phenotype_df[phenotype]
    # )
    res_dict = {}
    logger.info("Computing residual phenotypes")
    # for phenotype_suffix in ["", "_qt"]:     # comment this back in !!! 
    for phenotype_suffix in [""]:
        phenotype_col = f"{phenotype}{phenotype_suffix}"
        res_dict[phenotype_col] = {}

        for use_common_prs in [True, False]:
            covariate_cols = copy.deepcopy(base_covariate_cols)

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



@click.group()
def cli():
    pass


def find_substring_position(substring, string_list):
    matches = [index for index, string in enumerate(string_list) if substring in string]
    return matches[0] if matches else -1


@cli.command()
@click.option("--discoveries-file", type=click.Path(exists=True))
@click.option("--burdens-dir", type=click.Path(exists=True))
@click.option("--out-dir", type=click.Path(), default="./postprocessing/r_data")
@click.option("--debug", is_flag=True)
@click.option("--phenotype", type=str)
@click.option("--burden-type", multiple=True)  # , default = ('deeprvat', 'plof'))
@click.argument("sample_file_train", type=click.Path(exists=True))
@click.argument("sample_file_test", type=click.Path(exists=True))
@click.argument("config_file", type=click.Path(exists=True))
def prep_r_regression_data(
    discoveries_file,
    burdens_dir,
    out_dir,
    debug,
    phenotype,
    burden_type,
    sample_file_train,
    sample_file_test, 
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
    phenotype_df = phenotype_df.loc[[ i for i in phenotype_df.index if "W" not in i]]
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
    _, this_phenotype_df = add_corrected_phenotype_cols(phenotype, phenotype_df, prs_df, config)

    logger.info(this_phenotype_df.columns)

    discoveries = pd.read_parquet(discoveries_file)
    genes_this_pheno = compare_with_baseline_studies(
        phenotype,
        discoveries,
    )
    logger.info("Loading burdens")
    burdens_dict = {}
    for burden_type in burden_types:
        query_burden = burden_type.split("_")[0] # TODO can be deleted
        logger.info(f"Reading burdens for {query_burden}")

        with open(
            f"{burdens_dir}/{query_burden}_burdens.pickle", "rb"
        ) as f:
            this_burdens = pickle.load(f)
        burdens_dict[burden_type] = this_burdens
    

    sample_files = {'train': sample_file_train, 'test': sample_file_test}
    sample_dict = {}
    for split, file in sample_files.items():
        with open(file, 'rb') as f:
            this_samples = pickle.load(f)
            logger.info(f'Number of {split} samples: {len(this_samples)}')
            import ipdb; ipdb.set_trace()
            this_samples = list(set(this_samples).intersection(set(this_phenotype_df.index)).intersection(set(prs_df.dropna().index)))
            logger.info(f'Number of {split} samples after intersecting for samples with non-na pheno values: {len(this_samples)}')
            sample_dict[split] = this_samples

    for split in splits:
        data_dict = export_all_data(
            burdens_dict=burdens_dict,
            prs_df=prs_df,
            phenotype_df=this_phenotype_df,
            genes_this_pheno=genes_this_pheno,
            ordered_samples_to_keep=sample_dict[split],
            covariate_cols=config['covariates'],
            phenotype_cols=[
                col for col in this_phenotype_df.columns if phenotype in col
            ],
        )  
        logger.info("writing data")
        for btype, this_df in data_dict.items():
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
                            x_out_file = f"{out_dir}/{split}-{btype}-{gene_list}-{x_y_key}.parquet"
                            this_df.to_parquet(x_out_file)

        del data_dict

    logger.info("finished")


if __name__ == "__main__":
    cli()