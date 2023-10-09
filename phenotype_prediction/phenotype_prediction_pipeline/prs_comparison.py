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


PHENOTYPE_MAPPING = {
    "Forced_expiratory_volume_in_1_second_FEV1": "FEV1",
    "Mean_platelet_thrombocyte_volume": "MPTVS",
    "Body_mass_index_BMI": "BMI",
    "IGF_1": "IGF-1",
    "Red_blood_cell_erythrocyte_count": "Erythrocyte count",
}


def load_prs(prs_id, prs_file):
    logger.info(f"Retrieving prs from {prs_file} using PRS id {prs_id}")
    prs_df = pd.read_parquet(prs_file)[[prs_id]]
    prs_df.columns = ['common_PRS']
    return prs_df


def compare_with_baseline_studies(
    phenotype,
    discoveries,
    discoveries_all_samples=None,
):
    phenotype = phenotype.replace("_standardized", "")

    phenotype = (
        PHENOTYPE_MAPPING[phenotype]
        if phenotype in PHENOTYPE_MAPPING.keys()
        else " ".join(phenotype.split("_"))
    )
    logger.info(phenotype)
    genes_deeprvat_new = set(
        discoveries.query(
            'Trait == @phenotype & Method == "DeepRVAT" & \
    `Discovery type` == "New DeepRVAT discovery"'
        )["gene"]
    )
    genes_deeprvat_all = set(
        discoveries.query('Trait == @phenotype & Method == "DeepRVAT"')["gene"]
    )
    genes_baseline = set(
        discoveries.query('Trait == @phenotype  & Method == "Burden/SKAT combined"')[
            "gene"
        ]
    )
    if len(genes_baseline) == 0:
        top_n_genes = 1
        top_n_genes = min(len(genes_deeprvat_all), top_n_genes)
        logger.info(f'Using top {top_n_genes} of DeepRVAT discoveries for baseline because there were no baseline discoveries')
        genes_baseline = set(discoveries.query('Trait == @phenotype & Method == "DeepRVAT"')\
            .sort_values('pval')['gene'].drop_duplicates().head(top_n_genes))
    gene_dict = {
        "deeprvat_discoveries": genes_deeprvat_all,
        "deeprvat_novel": genes_deeprvat_new,
        "baseline_only": genes_baseline,
    }
    logger.info([f"{key} : {len(gene_dict[key])}" for key in gene_dict.keys()])

    if discoveries_all_samples is not None:
        genes_deeprvat_all_samples = set(
            discoveries_all_samples.query(
                'phenotype == @phenotype & Method == "DeepRVAT"'
            )["gene"]
        )
        gene_dict["deeprvat_all_samples"] = genes_deeprvat_all_samples
    logger.info(
        f"Number of baseline discoveries: {len(genes_baseline)} \n Number of DeepRVAT discoveries: {len(genes_deeprvat_all)} \n Number of Novel DeepRVAT discoveries: {len(genes_deeprvat_new)}"
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
        i for i in data_dict["deeprvat"]["baseline_only"][0]["x"].keys() if "gene_" in i
    ]
    corr_dict = {}
    for repeat in data_dict["deeprvat"]["baseline_only"].keys():
        corr_dict[f"repeat_{repeat}"] = {}
        for btype in data_dict["baseline"].keys():
            logger.info(f"Computing correlation for btype {btype}")
            corr_list = []
            for gene in genes:
                try:
                    scores_b = data_dict["baseline"][btype]["baseline_only"]["x"][gene]
                    scores_d = data_dict["deeprvat"]["baseline_only"][repeat]["x"][gene]
                    corr_list.append(round(stats.pearsonr(scores_d, scores_b)[0], 3))
                except:
                    logger.info(
                        f"Gene {gene} missing from deeprvat or baseline burden (likely {btype})"
                    )
            corr_dict[f"repeat_{repeat}"][btype] = corr_list
    avg_corr = {
        repeat: {
            btype: np.round(np.mean([np.abs(i) for i in corr_list]), 3)
            for btype, corr_list in subitem.items()
        }
        for repeat, subitem in corr_dict.items()
    }
    logger.info(f"Average correlation plof-deeprvat across all repeats: {avg_corr}")


def export_all_data(
    burdens_dict,
    prs_df,
    phenotype_df,
    genes_this_pheno,
    ordered_samples_to_keep,
    sample_mask,
    covariate_cols,
    phenotype_cols,
    baseline_sample_mask=None,
):
    data_dict = {"baseline": {}, "deeprvat": {}}

    # gene_lists = ['baseline_only', 'deeprvat_discoveries', 'genebass']
    gene_lists = ["baseline_only", "deeprvat_discoveries"]

    # 'deeprvat_all_samples']
    baseline_btypes = list(burdens_dict.keys())
    baseline_btypes.remove("deeprvat")
    data_dict["baseline"] = {
        btype: {gene_list: {} for gene_list in gene_lists} for btype in baseline_btypes
    }

    for btype in baseline_btypes:
        for gene_list in gene_lists:
            gene_list_suffix = "" if gene_list == "baseline_only" else "_deeprvat_genes"

            this_data_dict = get_model_data(
                burdens=burdens_dict[btype],
                btype=btype,
                genes_this_pheno=genes_this_pheno,
                gene_list=gene_list,
                sample_mask=sample_mask,
                ordered_sample_ids=ordered_samples_to_keep,
                prs_df=prs_df,
                phenotype_df=phenotype_df,
                phenotype_cols=phenotype_cols,
            )
            data_dict["baseline"][btype][gene_list] = this_data_dict

    for gene_list in gene_lists:
        this_data_dict = get_model_data(
            burdens=burdens_dict["deeprvat"],
            btype="deeprvat",
            genes_this_pheno=genes_this_pheno,
            gene_list=gene_list,
            sample_mask=sample_mask,
            ordered_sample_ids=ordered_samples_to_keep,
            prs_df=prs_df,
            phenotype_df=phenotype_df,
            phenotype_cols=phenotype_cols,
        )

        data_dict["deeprvat"][gene_list] = this_data_dict

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
    sample_mask,
    ordered_sample_ids,
    prs_df,
    phenotype_df,
    phenotype_cols,
):
    logger.info(f"Retreiveing model data for btype {btype}")
    query_genes = [
        gene for gene in genes_this_pheno[gene_list] if gene in list(burdens.keys())
    ]
    logger.info(
        f"Number of missing query genes: {len(query_genes) - len(genes_this_pheno[gene_list])}"
    )
    logger.info(
        f"Missing  genes: {set(genes_this_pheno[gene_list]) - set(burdens.keys())}"
    )

    if btype == "deeprvat":
        repeats = burdens[list(burdens.keys())[0]].keys()
        this_data_dict = {repeat: {} for repeat in repeats}
        # average predictions over all repeats
        for repeat in repeats:
            this_deeprvat_burdens = [burdens[gene][repeat] for gene in query_genes]
            this_deeprvat_burdens = np.column_stack(this_deeprvat_burdens)[sample_mask]

            logger.info(f"deeprvat burden shape {this_deeprvat_burdens.shape}")

            y, x = _get_model_data(
                gene_burdens=this_deeprvat_burdens,
                ordered_sample_ids=ordered_sample_ids,
                this_prs=prs_df,
                phenotype_df=phenotype_df,
                gene_ids=query_genes,
                phenotype_cols=phenotype_cols,
            )

            this_data_dict[repeat] = {"y": y, "x": x}
    else:
        this_baseline_burdens = [burdens[gene] for gene in query_genes]
        this_baseline_burdens = np.column_stack(this_baseline_burdens)[sample_mask]
        logger.info(f"Baseline burdens shape: {this_baseline_burdens.shape}")
        y, x = _get_model_data(
            gene_burdens=this_baseline_burdens,
            ordered_sample_ids=ordered_sample_ids,
            this_prs=prs_df,
            phenotype_df=phenotype_df,
            gene_ids=query_genes,
            phenotype_cols=phenotype_cols,
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
    # import ipdb; ipdb.set_trace()
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
