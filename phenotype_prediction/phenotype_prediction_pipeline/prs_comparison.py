import pandas as pd
from pathlib import Path
import zarr
import numpy as np
from scipy.sparse import coo_matrix, spmatrix
from scipy import stats
import statsmodels.api as sm
import pickle
from deeprvat.data import DenseGTDataset


import sys
import os
import logging

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

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

PHENOTYPE_MAPPING = {
    "Forced_expiratory_volume_in_1_second_FEV1": "FEV1",
    "Mean_platelet_thrombocyte_volume": "MPTVS",
    "Body_mass_index_BMI": "BMI",
    "IGF_1": "IGF-1",
    "LDL_direct_statin_corrected": 'LDL direct',
    'Cholesterol_statin_corrected': 'Cholesterol',
    "Red_blood_cell_erythrocyte_count": "Erythrocyte count",
}


def load_prs(prs_id, prs_file):
    logger.info(f"Retrieving prs from {prs_file} using PRS id {prs_id}")
    prs_df = pd.read_parquet(prs_file)[[prs_id]]
    prs_df.columns = ["common_PRS"]
    return prs_df


def compare_with_baseline_studies(
    phenotype,
    discoveries,
):
    logger.info(phenotype)

    genes_deeprvat_all = set(
        discoveries.query('phenotype == @phenotype & Method == "DeepRVAT"')["gene"]
    )
    genes_baseline = set(
        discoveries.query('phenotype == @phenotype  & Method == "Burden/SKAT combined"')[
            "gene"
        ]
    )

    gene_dict = {
        "deeprvat_discoveries": genes_deeprvat_all,
        "baseline_only": genes_baseline,
    }
    logger.info([f"{key} : {len(gene_dict[key])}" for key in gene_dict.keys()])


    logger.info(
        f"Number of baseline discoveries: {len(genes_baseline)} \n Number of DeepRVAT discoveries: {len(genes_deeprvat_all)}"
    )

    return gene_dict


def get_vars_pers_sample_and_gene(gene, G_full, grouped_annotations):
    try:
        annotation_df = grouped_annotations.get_group(gene)
        variant_ids = annotation_df.index.unique().to_numpy()
        G = G_full.tocsc()[:, variant_ids].todense()
        G = np.asarray(G)
        vars_per_sample = np.sum(G, axis=1)
        samples_with_variant = vars_per_sample[vars_per_sample > 0].shape[0]
        # logger.info(f'Number of samples with a variant: {samples_with_variant}')
    except:
        # do this if grouped_annotations.get_group(gene) returns an error since there is no plof variant in the gene
        logger.info(f"No variants in this variant category in the gene")
        vars_per_sample = np.zeros((G_full.shape[0],))
    return vars_per_sample


def _get_burden_correlation(data_dict):
    genes = [
        i for i in data_dict["deeprvat"]["baseline_only"]["x"].keys() if "gene_" in i
    ]
    corr_dict = {}
    for btype in data_dict.keys():
        logger.info(f"Computing correlation for btype {btype}")
        corr_list = []
        for gene in genes:
            try:
                scores_b = data_dict[btype]["baseline_only"]["x"][gene]
                scores_d = data_dict["deeprvat"]["baseline_only"]["x"][gene]
                corr_list.append(round(stats.pearsonr(scores_d, scores_b)[0], 3))
            except:
                logger.info(
                    f"Gene {gene} missing from deeprvat or baseline burden (likely {btype})"
                )
        corr_dict[btype] = corr_list
    avg_corr = {
            btype: np.round(np.mean([np.abs(i) for i in corr_list]), 3)
            for btype, corr_list in corr_dict.items()
    }
    logger.info(f"Average correlation baseline-deeprvat across all repeats: {avg_corr}")


def export_all_data(
    burdens_dict,
    prs_df,
    phenotype_df,
    genes_this_pheno,
    ordered_samples_to_keep,
    covariate_cols,
    phenotype_cols,
):

    gene_lists = ["baseline_only", "deeprvat_discoveries"]

    data_dict = {
        btype: {gene_list: {} for gene_list in gene_lists} for btype in burdens_dict.keys()
    }
    for btype in burdens_dict.keys():
        for gene_list in gene_lists:

            this_data_dict = get_model_data(
                burdens=burdens_dict[btype],
                btype=btype,
                genes_this_pheno=genes_this_pheno,
                gene_list=gene_list,
                ordered_sample_ids=ordered_samples_to_keep,
                prs_df=prs_df,
                phenotype_df=phenotype_df,
                phenotype_cols=phenotype_cols,
                covariate_cols=covariate_cols
            )
            data_dict[btype][gene_list] = this_data_dict       
    ########################################################
    # sanity chek: correlation of plof scores with deeprvat
    logger.info("Computing correlation of different burden types for significant genes")
    _get_burden_correlation(data_dict)

    return data_dict


def get_model_data(
    burdens,
    btype,
    genes_this_pheno,
    gene_list,
    ordered_sample_ids,
    prs_df,
    phenotype_df,
    phenotype_cols,
    covariate_cols
):
    logger.info(f"Retreiveing model data for btype {btype}")
    query_genes = [
        gene for gene in genes_this_pheno[gene_list] if gene in list(burdens.keys())
    ]
    logger.info(
        f"Number of missing query genes: {len(genes_this_pheno[gene_list])- len(query_genes)}"
    )
    logger.info(
        f"Missing  genes: {set(genes_this_pheno[gene_list]) - set(burdens.keys())}"
    )

    this_burdens = [burdens[gene] for gene in query_genes]
    burden_sample_map, _ = get_matched_sample_indices(burdens['sample_ids'], ordered_sample_ids)
    # mapping sanity check, can be removed in future 
    for i in np.random.randint(len(ordered_sample_ids), size=10):

        assert burdens['sample_ids'][burden_sample_map[i]] == ordered_sample_ids[i]
    test = np.column_stack(this_burdens)  # sanity check, can be removed in future 
    this_burdens = np.column_stack(this_burdens)[burden_sample_map]
    # sanity check, can be removed in future 

    for  i in np.random.choice(np.where(this_burdens.sum(axis =1) > 0)[0], 10):
        assert np.all(test[burden_sample_map[i], :] == this_burdens[i, :])

    logger.info(f"Baseline burdens shape: {this_burdens.shape}")
    y, x = _get_model_data(
        gene_burdens=this_burdens,
        ordered_sample_ids=ordered_sample_ids,
        this_prs=prs_df,
        phenotype_df=phenotype_df,
        gene_ids=query_genes,
        phenotype_cols=phenotype_cols,
        covariate_cols=covariate_cols,

    )
    this_data_dict = {"y": y, "x": x}

    return this_data_dict


def _get_model_data(
    gene_burdens,
    ordered_sample_ids,
    this_prs,
    phenotype_df,
    phenotype_cols,
    covariate_cols=["age", "genetic_sex", *[f"genetic_PC_{i}" for i in range(1, 21)]],
    gene_ids=None,
):

    if gene_burdens is not None:
        X = np.column_stack(
            [
                this_prs.loc[ordered_sample_ids],
                gene_burdens,
                phenotype_df.loc[ordered_sample_ids][covariate_cols],
            ]
        )

        if gene_ids is None:
            gene_cols = [f"gene_{i}" for i in range(gene_burdens.shape[1])]
        else:
            gene_cols = [f"gene_{gene_id}" for gene_id in gene_ids]
        colnames = ["prs"] + gene_cols + covariate_cols

        # *[f'gene_{i}' for i in range(gene_burdens.shape[1])], *covariate_cols]

    else:
        X = np.column_stack(
            [
                this_prs.loc[ordered_sample_ids],
                phenotype_df.loc[ordered_sample_ids][covariate_cols],
            ]
        )
        colnames = ["prs", *covariate_cols]
    X = pd.DataFrame(X, columns=colnames, index=ordered_sample_ids)
    Y = phenotype_df.loc[ordered_sample_ids][phenotype_cols]

    return Y, X


def get_residual_phenotypes(
    gene_burdens,
    ordered_sample_ids,
    phenotype_df,
    covariate_cols,
    phenotype_col,
    substract_covariates=True,
    fitted_model=None,
):
    if gene_burdens is not None:
        # X = np.column_stack([this_prs.loc[ordered_sample_ids],
        #                     gene_burdens,
        #                     phenotype_df.loc[ordered_sample_ids][covariate_cols]])
        # colnames = ['prs', *[f'gene_{i}' for i in range(gene_burdens.shape[1])],
        #             *covariate_cols]
        X = np.column_stack(
            [gene_burdens, phenotype_df.loc[ordered_sample_ids][covariate_cols]]
        )
        colnames = [
            *[f"gene_{i}" for i in range(gene_burdens.shape[1])],
            *covariate_cols,
        ]
    else:
        X = np.column_stack([phenotype_df.loc[ordered_sample_ids][covariate_cols]])
        colnames = covariate_cols
    # import ipdb; ipdb.set_trace()
    X = pd.DataFrame(X, columns=colnames, index=ordered_sample_ids)
    Y = phenotype_df.loc[ordered_sample_ids][phenotype_col]
    X = sm.add_constant(X)
    # TODO: compute correlations here
    if fitted_model is None:
        logger.info("Fitting regression model")
        model = sm.OLS(Y, X)
        results = model.fit()
    else:
        logger.info("Computing PRS based on provided model parameters")
        results = fitted_model
    betas = results.params
    Y_hat = X * betas
    Y_hat = Y_hat.sum(axis=1)

    residuals = Y - Y_hat

    return Y_hat, residuals, results
