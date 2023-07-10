import pickle
import math
import sys
import os
import numpy as np
import itertools
import logging
import pandas as pd
import subprocess
from pathlib import Path
from statsmodels.stats.multitest import fdrcorrection
import yaml

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


PVAL_CORRECTION_TYPES = {
    "FDR": "",
    "Bonferroni": "_bf",
}

config_cols = [
    "prop_causal_variants",
    "max_rare_af",
    "var_genetic_effects",
    "var_noise",
    "var_expl_variant_noise",
    "var_expl_cont",
    "var_expl_binary",
    "var_expl_maf",
    "factor_var_expl_cont",
]

column_renamer_reverse = {
    "vprop": "prop_causal_variants",
    "maxaf": "max_rare_af",
    "vnoise": "var_noise",
    "vgen": "var_genetic_effects",
    "vvarnoise": "var_expl_variant_noise",
    "vcont": "var_expl_cont",
    "vbin": "var_expl_binary",
    "vmaf": "var_expl_maf",
    "vextra": "var_extra_group",
}

BASELINE_METHODS_DICT = {
    "Burden pLOF": ["plof-burden"],
    "Burden missense": ["missense-burden"],
    "SKAT pLOF": ["plof-skat"],
    "SKAT missense": ["missense-skat"],
    "Burden/SKAT combined": [
        "missense-burden",
        "missense-skat",
        "plof-burden",
        "plof-skat",
    ],
}


def collect_all_results(
    exp_name,
    quantiles=[None],
    caf_threshold=None,
    local_res_dir="./",
    config_cols=config_cols,
    alpha=0.05,
    return_associations=False,
    genes_to_exclude=None,
):
    local_base_target_dir = Path(f"{local_res_dir}/{exp_name}")

    simulation_config = pd.read_parquet(
        f"{local_base_target_dir}/simulation_config.parquet"
    )
    if "n_total_causal_variants" in simulation_config.columns:
        simulation_config["n_total_causal_variants"] = simulation_config[
            "n_total_causal_variants"
        ].astype(int)
    simulation_config_complete = pd.read_parquet(
        f"{local_base_target_dir}/simulation_config_complete.parquet"
    )
    if "n_total_causal_variants" in simulation_config_complete.columns:
        simulation_config_complete[
            "n_total_causal_variants"
        ] = simulation_config_complete["n_total_causal_variants"].astype(int)

    n_configs = len(simulation_config)
    exp_configs = [
        [f"{key}:{simulation_config[key].iloc[i]}" for key in simulation_config.keys()]
        for i in range(n_configs)
    ]
    config_folders = ["_".join(sublist) for sublist in exp_configs]
    config_names = [";".join(sublist) for sublist in exp_configs]
    res_list = []
    association_list = []
    for i in range(len(exp_configs)):
        for q in quantiles:
            res_df, corrected_associations = get_all_metrics(
                local_base_target_dir,
                config_folders[i],
                q,
                caf_threshold,
                alpha,
                genes_to_exclude,
            )
            res_df = (
                res_df.assign(**simulation_config.iloc[i].to_dict())
                .assign(config=config_names[i])
                .merge(simulation_config_complete, how="left")
                .assign(quantile=q)
            )
            corrected_associations = (
                corrected_associations.assign(**simulation_config.iloc[i].to_dict())
                .assign(config=config_names[i])
                .merge(simulation_config_complete, how="left")
                .assign(quantile=q)
            )

            res_list.append(res_df)
            association_list.append(corrected_associations)

    all_corrected_associations = pd.concat(association_list)

    all_res_df = (
        pd.concat(res_list)
        .assign(
            sig_line_alpha=lambda df: df["metric"].map(
                lambda m: 1.0 if m in ["FDR", "Power"] else 0.0
            )
        )
        .assign(
            sig_line_value=lambda df: df["metric"].map(
                lambda m: alpha if m in ["FDR", "FDR_bf"] else 0
            )
        )
    )

    categories = list(BASELINE_METHODS_DICT.keys()) + [
        "DeepRVAT_wo_baseline",
        "DeepRVAT",
    ]
    test_type_cat = pd.Categorical(all_res_df["test"], categories=categories)
    all_res_df = all_res_df.assign(test=test_type_cat)

    avg_metrics = (
        all_res_df[["metric", "test", "config", "quantile", "value"]]
        .groupby(["metric", "test", "config", "quantile"])
        .agg({"value": [np.mean, np.std]})["value"]
        .reset_index()
        .merge(all_res_df[config_cols + ["config"]].drop_duplicates())
    )
    avg_metrics_display = (
        avg_metrics[["metric", "test", "config", "quantile", "mean", "std"]]
        .query('metric in ["Power", "FDR", "TP", "FP"]')
        .sort_values(["config", "metric", "test"])
    )

    if "vmaf" in all_res_df.columns:
        all_res_df = all_res_df.rename(columns=column_renamer_reverse)
        avg_metrics = avg_metrics.rename(columns=column_renamer_reverse)

    if return_associations:
        return all_res_df, avg_metrics, all_corrected_associations
    else:
        return all_res_df, avg_metrics


def get_all_metrics(
    exp_dir,
    exp_config,
    quantile=None,
    caf_threshold=None,
    alpha=0.05,
    genes_to_exclude=None,
):
    try:
        all_associations = pd.read_parquet(
            f"{exp_dir}/{exp_config}/all_associations.parquet"
        )
    except:
        all_associations = pd.read_parquet(
            f"{exp_dir}/{exp_config}/eval/all_associations.parquet"
        )
    sim_repeats = [i for i in all_associations["repeat"].unique() if "sim_repeat_" in i]
    print(sim_repeats)
    all_metrics = []
    all_corrected_associations = []
    for sim_repeat in sim_repeats:
        associations_this_sim_repeat = all_associations.copy().query(
            "repeat == @sim_repeat"
        )
        this_metrics, this_corrected_results = get_metrics_for_sim_repeat(
            associations_this_sim_repeat,
            quantile,
            caf_threshold=caf_threshold,
            alpha=alpha,
            correlated_genes_to_remove_dict=genes_to_exclude,
        )
        all_metrics.append(pd.DataFrame(this_metrics))
        all_corrected_associations.append(this_corrected_results)
    metrics_df = (
        pd.concat(all_metrics, keys=sim_repeats)
        .reset_index()
        .rename(columns={"level_0": "repeat", "level_1": "metric"})
    )

    all_corrected_associations_df = pd.concat(all_corrected_associations)
    metrics_df = metrics_df.melt(
        id_vars=["repeat", "metric"], value_name="value", var_name="test"
    )
    return metrics_df, all_corrected_associations_df


def get_metrics_for_sim_repeat(
    associations_this_sim_repeat,
    quantile=None,
    caf_threshold=None,
    n_deeprvat_repeats=6,
    n_causal_genes=100,
    min_sig_deeprvat_repeats=1,
    alpha=0.05,
    correlated_genes_to_remove_dict=None,
):
    ########
    metrics = {}

    baseline_methods = BASELINE_METHODS_DICT["Burden/SKAT combined"]
    # baseline results are already EAC filtered

    baseline_res = associations_this_sim_repeat.copy().query(
        'method in @baseline_methods & correction_method == "FDR"'
    )
    baseline_res[["vtype", "ttype"]] = baseline_res["method"].str.split(
        "-", 1, expand=True
    )
    if correlated_genes_to_remove_dict == None:
        logger.info("not excluding correlated genes")
        baseline_res["keep_gene"] = True
    else:
        correlated_genes_to_remove = list(
            correlated_genes_to_remove_dict[quantile]["all_vtypes"]
        )
        # logger.info(f'Number of genes to remove {len(correlated_genes_to_remove)}')
        baseline_res = baseline_res.query(
            "gene not in  @correlated_genes_to_remove"
        ).assign(keep_gene=True)
    logger.info(
        f'Total number of genes kept in baseline: {len(baseline_res["gene"].unique())}'
    )

    baseline_res_corrected = pval_correct_all_methods(
        baseline_res, PVAL_CORRECTION_TYPES, alpha=alpha
    )

    for test_name, tests in BASELINE_METHODS_DICT.items():
        metrics[test_name] = get_metrics_(
            baseline_res_corrected.query("method in @tests"), n_causal_genes
        )

    deeprvat_res = associations_this_sim_repeat.query(
        f'experiment == "DeepRVAT ({n_deeprvat_repeats} repeats)"'
    ).query('method == "DeepRVAT" & correction_method == "FDR"')
    if quantile is not None:
        # also remove correlated genes from DeepRVAT to use the same gene set everywhere
        deeprvat_res = deeprvat_res.query("gene not in  @correlated_genes_to_remove")
    deeprvat_res_corrected = pval_correct_all_methods(
        deeprvat_res, PVAL_CORRECTION_TYPES, alpha=alpha
    )
    # Filter for genes significant in >1 DeepRVAT repeat
    genes_to_keep_deeprvat = (
        deeprvat_res_corrected.query('correction_method == "FDR"')[
            ["gene", "significant"]
        ]
        .groupby("gene")
        .sum()
        .query("significant > @min_sig_deeprvat_repeats")
        .index.to_list()
    )
    deeprvat_res_corrected = deeprvat_res_corrected.assign(
        keep_gene=lambda df: df["gene"].map(
            lambda gene: True if gene in genes_to_keep_deeprvat else False
        )
    )
    deeprvat_res = deeprvat_res.assign(
        keep_gene=lambda df: df["gene"].map(
            lambda gene: True if gene in genes_to_keep_deeprvat else False
        )
    )
    metrics["DeepRVAT_wo_baseline"] = get_metrics_(
        deeprvat_res_corrected.query("keep_gene"), n_causal_genes
    )

    baseline_deeprvat_combined = pd.concat(
        [deeprvat_res, baseline_res.query("keep_gene")]
    )
    baseline_deeprvat_combined_corrected = pval_correct_all_methods(
        baseline_deeprvat_combined, PVAL_CORRECTION_TYPES, alpha=alpha
    )

    metrics["DeepRVAT"] = get_metrics_(
        baseline_deeprvat_combined_corrected.query("keep_gene"), n_causal_genes
    )

    all_corrected_results = (
        pd.concat(
            [
                baseline_res_corrected,
                deeprvat_res_corrected,
                baseline_deeprvat_combined_corrected,
            ],
            keys=["Baseline", "DeepRvat_wo_baseline", "DeepRVAT"],
        )
        .reset_index(level=0)
        .rename(columns={"level_0": "combination"})
    )

    return metrics, all_corrected_results


def get_metrics_(this_associations, n_causal=100):
    this_metric = {}

    for correction_type, sig_col_suffix in PVAL_CORRECTION_TYPES.items():
        # logger.info(f'Correction type: {correction_type}')
        tp = this_associations.query(
            f"correction_method == @correction_type & significant & causal_gene"
        ).drop_duplicates("gene")
        # drop duplicates to account for the fact the genes might pop up multiple types
        # if significant in multiple baseline tests/repeats

        n_tp = len(tp)
        this_metric[f"TP{sig_col_suffix}"] = n_tp
        # logger.info(f'Number of true positives: {len(tp)}')

        power = n_tp / n_causal
        this_metric[f"Power{sig_col_suffix}"] = power
        # logger.info(f'Power all causal genes: {power:.4f}')

        fp = this_associations.query(
            f"correction_method == @correction_type &  significant & causal_gene == False"
        ).drop_duplicates("gene")
        n_fp = len(fp)
        this_metric[f"FP{sig_col_suffix}"] = n_fp
        # logger.info(f'Number of false-positives: {n_fp}')
        if n_fp > 0 or n_tp > 0:
            fdr = n_fp / (n_fp + n_tp)

            this_metric[f"FDR{sig_col_suffix}"] = fdr
            # logger.info(f'FDR: {fdr:.4f}')
        else:
            this_metric[f"FDR{sig_col_suffix}"] = 0

        # logger.info('True-positives:')
        # print(tp)

    return this_metric


def fdrcorrect_df(group: pd.DataFrame, alpha: float) -> pd.DataFrame:
    group = group.copy()

    rejected, pval_corrected = fdrcorrection(group["pval"], alpha=alpha)
    group["significant"] = rejected
    group["pval_corrected"] = pval_corrected
    return group


def bfcorrect_df(group: pd.DataFrame, alpha: float) -> pd.DataFrame:
    group = group.copy()

    pval_corrected = group["pval"] * len(group)
    group["significant"] = pval_corrected < alpha
    group["pval_corrected"] = pval_corrected
    return group


def pval_correction(group: pd.DataFrame, alpha: float, correction_type: str = "FDR"):
    if correction_type == "FDR":
        return fdrcorrect_df(group, alpha)
    elif correction_type == "Bonferroni":
        return bfcorrect_df(group, alpha)
    # elif correction_type == 'CCT/FDR':
    #     return cctcorrect_df(group, alpha)
    # elif correction_type == 'CCT/Bonferroni':
    #     return cctbfcorrect_df(group, alpha)
    else:
        raise ValueError(
            f"Unknown correction type: {correction_type}. "
            "Valid values are '(CCT/)FDR' and '(CCT/)Bonferroni'."
        )


def pval_correct_all_methods(result, pval_correction_types, alpha=0.05):
    result = result.loc[:, ~result.columns.str.endswith("_bf")]

    corrected_results = []
    for correction_type in pval_correction_types.keys():
        this_corrected_result = pval_correction(
            result.copy(), alpha, correction_type=correction_type
        ).assign(correction_method=correction_type)

        this_corrected_result["-log10pval_corrected"] = -np.log10(
            this_corrected_result["pval_corrected"]
        )

        corrected_results.append(this_corrected_result)

    return pd.concat(corrected_results)


################### additional functions ############################################################
##################################################################################################################
def get_maf_filter_var_explained(experiment_dir):
    logger.info("Reading variant properties from {experiment_dir}")

    exp_renamer = {
        "vprop:0.2_maxaf:0.01_vgen:0.1_vnoise:0.95_vvarnoise:0.1_vcont:0.3_vbin:0.3_vmaf:0.3_expbin:1": "0.01%",
        "vprop:0.2_maxaf:0.01_vgen:0.1_vnoise:0.95_vvarnoise:0.1_vcont:0.4_vbin:0.4_vmaf:0.1_expbin:10": "1%",
        "vprop:0.2_maxaf:0.01_vgen:0.1_vnoise:0.95_vvarnoise:0.1_vcont:0.35_vbin:0.35_vmaf:0.2_expbin:3": "0.1%",
    }

    config_file = f"{experiment_dir}/config.yaml"
    with open(config_file) as f:
        config = yaml.safe_load(f)

    logger.info("Reading annotation data frame")

    annotation_file = config["simulation"]["data"]["dataset_config"]["annotation_file"]
    columns = config["simulation"]["data"]["dataset_config"]["annotations"] + ["id"]

    annotation_df = pd.read_parquet(annotation_file, columns=columns).set_index("id")
    annotation_df

    n_sim_repeats = 10
    all_configs = os.listdir(f"{experiment_dir}")
    all_configs = [config for config in all_configs if "vprop" in config]

    logger.info("Collecting data")
    var_scores_dict = {config: {} for config in all_configs}
    for config in all_configs:
        for sim_repeat in range(n_sim_repeats):
            this_data_dir = (
                f"{experiment_dir}/{config}/sim_repeat_{sim_repeat}/simulation_data/"
            )
            this_var_scores = (
                pd.read_parquet(f"{this_data_dir}/variant_scores.parquet")
                .query("sim_causal")
                .set_index("id")["variant_score"]
            )
            var_scores_dict[config][sim_repeat] = this_var_scores

    stacked_bins = [0, 0.0001, 0.001, 0.01]
    bins = stacked_bins

    n_samples_ukbb_cohort = 167245
    var_weights_per_bin_dict = {}

    logger.info("Computing variance explained by each MAF bin ")

    weight_df_list = []
    for weight_fnc in var_scores_dict.keys():
        for sim_repeat, this_df in var_scores_dict[weight_fnc].items():
            this_df = var_scores_dict[weight_fnc][sim_repeat].copy()
            this_df = (
                this_df.to_frame()
                .join(annotation_df.loc[this_df.index]["UKB_AF"], how="left")
                .assign(weight_fnc=weight_fnc)
                .assign(sim_repeat=sim_repeat)
            )
            weight_df_list.append(this_df)

    weighted_vars_df = pd.concat(weight_df_list)
    weighted_vars_df["allele_count"] = (
        weighted_vars_df["UKB_AF"] * 2 * n_samples_ukbb_cohort
    )
    weighted_vars_df["effect_size"] = (
        weighted_vars_df["allele_count"] * weighted_vars_df["variant_score"]
    )

    summed_effect_size = (
        weighted_vars_df.groupby(["weight_fnc", "sim_repeat"])["effect_size"]
        .sum()
        .reset_index()
        .rename(columns={"effect_size": "summed_effect_size"})
    )
    var_weights_per_bin = (
        weighted_vars_df.groupby(
            [pd.cut(weighted_vars_df.UKB_AF, bins)] + ["weight_fnc", "sim_repeat"]
        )["effect_size"]
        .sum()
        .reset_index()
        .rename(columns={"UKB_AF": "interval"})
        .merge(summed_effect_size, how="left")
    )
    var_weights_per_bin["rel_effect"] = (
        var_weights_per_bin["effect_size"] / var_weights_per_bin["summed_effect_size"]
    )
    var_weights_per_bin = (
        var_weights_per_bin.groupby(["interval", "weight_fnc"])
        .mean()
        .reset_index()
        .drop(columns="sim_repeat")
        .assign(
            group=lambda df: df["weight_fnc"].map(lambda config: exp_renamer[config])
        )
    )
    var_weights_per_bin["interval"] = var_weights_per_bin["interval"].astype(str)

    logger.info("Saving data")

    return var_weights_per_bin
