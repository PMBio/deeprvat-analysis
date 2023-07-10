import pandas as pd
import numpy as np
import yaml
import sys
import logging
import click
import pickle
import itertools
import re
from plotnine import *
import copy
from pathlib import Path
from statsmodels.stats.multitest import fdrcorrection

from deeprvat.utils import pval_correction


logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

PVAL_CORRECTION_TYPES = {
    "FDR": "",
    "Bonferroni": "_bf",
    "CCT/FDR": "_cct",
    "CCT/Bonferroni": "_cct_bf",
}

correction_methods = ["FDR", "Bonferroni"]

pval_correction_types = {
    k: v for k, v in PVAL_CORRECTION_TYPES.items() if k in correction_methods
}

cols_to_keep = [
    "gene",
    "pval",
    "phenotype",
    "-log10pval",
    "significant",
    "pval_corrected",
    "seed_gene",
    "causal_gene",
    "method",
    "correction_method",
]

#################################################################################################
####### CAUTION!!!!!!
####### CAUTION!!!!!!
# The metrics ( metrics.pkl) computed here are not correct!!!!! ####
# Only use the aggregated associations (all_associations.parquet) returned by this script
# and then continue with the post-processing using the notebook provided by Eva!!!!
#################################################################################################


def pval_correct_all_methods(result, pval_correction_types, alpha=0.05):
    print("doing pval correction")

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


def add_simulated_causal_gene_col_(
    this_associations, sim_causal_gene_ids, seed_gene_ids=None
):
    this_associations = this_associations.assign(
        causal_gene=lambda df: df["gene"].map(
            lambda this_id: True if this_id in sim_causal_gene_ids else False
        )
    )
    if seed_gene_ids is not None:
        this_associations = this_associations.assign(
            seed_gene=lambda df: df["gene"].map(
                lambda this_id: True if this_id in seed_gene_ids else False
            )
        ).drop(columns="in_seed_genes")
    else:
        this_associations["seed_gene"] = False

    return this_associations


def get_metrics_(this_associations, n_causal, n_seed):
    pval_correction_types = {"FDR": "", "Bonferroni": "_bf"}  # ,
    this_metric = {}

    n_caual_non_seed = n_causal - n_seed
    for correction_type, sig_col_suffix in pval_correction_types.items():
        logger.info(f"Correction type: {correction_type}")
        tp = this_associations.query(
            f"correction_method == @correction_type & significant & causal_gene"
        ).drop_duplicates("gene")
        # drop duplicates to account for the fact the genes might pop up multiple types
        # if significant in multiple baseline tests/repeats

        n_tp = len(tp)
        this_metric[f"TP{sig_col_suffix}"] = n_tp
        logger.info(f"Number of true positives: {len(tp)}")

        tp_seed = this_associations.query(
            f"correction_method == @correction_type & significant & seed_gene"
        ).drop_duplicates("gene")
        this_metric[f"TP-seed{sig_col_suffix}"] = len(tp_seed)

        power = n_tp / n_causal
        this_metric[f"power{sig_col_suffix}"] = power
        logger.info(f"Power all causal genes: {power:.4f}")

        fp = this_associations.query(
            f"correction_method == @correction_type &  significant & ~causal_gene"
        ).drop_duplicates("gene")
        n_fp = len(fp)
        this_metric[f"FP{sig_col_suffix}"] = n_fp
        logger.info(f"Number of false-positives: {n_fp}")
        if n_fp > 0 or n_tp > 0:
            fdr = n_fp / (n_fp + n_tp)

            this_metric[f"FDR{sig_col_suffix}"] = fdr
            logger.info(f"FDR: {fdr:.4f}")
        else:
            this_metric[f"FDR{sig_col_suffix}"] = 0

        logger.info("True-positives:")
        print(tp)

    return this_metric


def compute_metrics_average_(metrics):
    repeats = copy.deepcopy(list(metrics.keys()))
    repeats = [repeat for repeat in repeats if "average" not in repeat]
    metrics["average"] = {}
    metrics["sd"] = {}
    methods = list(metrics[repeats[0]].keys())
    for method in methods:
        metrics["average"][method] = {}
        metrics["sd"][method] = {}
        eval_fields = list(metrics[repeats[0]][method].keys())
        for eval_field in eval_fields:
            metrics["average"][method][eval_field] = {}
            metrics["sd"][method][eval_field] = {}
            metric_types = list(metrics[repeats[0]][method][eval_field].keys())
            for metric_type in metric_types:
                values = [
                    metrics[repeat][method][eval_field][metric_type]
                    for repeat in repeats
                ]
                metrics["average"][method][eval_field][metric_type] = np.mean(values)
                metrics["sd"][method][eval_field][metric_type] = np.std(values)

    return metrics


def evaluate_baseline_(
    exp_name, repeat, baseline_test_combis, sim_causal_gene_ids, pval_correction_types
):
    res_base_dir = f"{exp_name}/sim_repeat_{repeat}/baseline/sim_phenotype/"

    associations = pd.read_parquet(
        f"{res_base_dir}/eval/burden_associations_testing.parquet"
    ).query("CAF > 50")

    associations = pval_correct_all_methods(associations, pval_correction_types)
    associations = add_simulated_causal_gene_col_(
        associations, sim_causal_gene_ids, None
    )[cols_to_keep]

    metrics = {}

    for test in associations["method"].unique():
        logger.info(f"Getting metrics for test {test}")
        metrics[test] = get_metrics_(
            associations.query("method == @test"),
            n_causal=len(sim_causal_gene_ids),
            n_seed=0,
        )
    metrics["baseline_union"] = get_metrics_(
        associations, n_causal=len(sim_causal_gene_ids), n_seed=0
    )

    return metrics, associations


def evaluate_deeprvat_(exp_name, repeat, sim_causal_gene_ids):
    test_name = "deeprvat"
    logger.info(f"analyzing resuts for {test_name}")
    result_dir = f"{exp_name}/sim_repeat_{repeat}/deeprvat/sim_phenotype/deeprvat/eval"

    seed_gene_ids = pd.read_parquet(f"{result_dir}/seed_genes.parquet")["id"].to_list()
    associations = (
        pd.read_parquet(f"{result_dir}/burden_associations_testing.parquet")
        .drop(columns=["in_baseline"])
        .query(
            'experiment == "DeepRVAT (6 repeats)" & \
           method == "DeepRVAT"'
        )
    )
    associations = add_simulated_causal_gene_col_(
        associations, sim_causal_gene_ids, seed_gene_ids
    )

    metrics = {}
    genes_to_keep = (
        associations.query('correction_method == "FDR"')[["gene", "significant"]]
        .groupby("gene")
        .sum()
        .query("significant > 1")
        .index.to_list()
    )
    metrics["DeepRVAT"] = get_metrics_(
        associations.query("gene in @genes_to_keep"),
        n_causal=len(sim_causal_gene_ids),
        n_seed=len(seed_gene_ids),
    )
    # genes_to_keep_with_baseline = associations.query('method == "Seed gene discovery" & significant \
    #     & experiment == "DeepRVAT (6 repeats)" & correction_method == "FDR"')['gene'].unique()
    # genes_to_keep_with_baseline = list(set(genes_to_keep_with_baseline).union(set(genes_to_keep)))
    # metrics['DeepRVAT_with_baseline'] = get_metrics_(associations.query('gene in @genes_to_keep_with_baseline & \
    #     experiment == "DeepRVAT (6 repeats)"'),
    #     n_causal=len(sim_causal_gene_ids), n_seed=len(seed_gene_ids))

    return metrics, associations


def evaluate_(exp_name, n_repeats, baseline_test_combis):
    all_associations = []
    all_metrics = {}
    for repeat in range(n_repeats):
        logger.info(f"Analyzing results for repeat {repeat}")
        sim_causal_gene_file = f"{exp_name}/sim_repeat_{repeat}/simulation_data/causal_genes_simulated.parquet"
        sim_causal_genes = pd.read_parquet(sim_causal_gene_file)
        sim_causal_gene_ids = sim_causal_genes["id"].to_list()

        all_metrics[f"sim_repeat_{repeat}"] = {}
        b_metrics, b_associations = evaluate_baseline_(
            exp_name,
            repeat,
            baseline_test_combis,
            sim_causal_gene_ids,
            pval_correction_types,
        )
        d_metrics, d_associations = evaluate_deeprvat_(
            exp_name, repeat, sim_causal_gene_ids
        )

        all_metrics[f"sim_repeat_{repeat}"]["baseline"] = b_metrics
        all_metrics[f"sim_repeat_{repeat}"]["deeprvat"] = d_metrics

        associations = pd.concat(
            [
                b_associations.assign(repeat=f"sim_repeat_{repeat}"),
                d_associations.assign(repeat=f"sim_repeat_{repeat}"),
            ]
        )

        all_associations.append(associations)

    all_metrics = compute_metrics_average_(all_metrics)
    all_associations = pd.concat(all_associations)

    return all_metrics, all_associations


##########
######### Evaluate the simulated data set


def _get_causal_variant_properties(
    exp_name, n_repeats, annotation_df, vars_to_ensgid_mapping
):
    causal_vars_with_anno_list = []
    causal_vars_per_gene_list = []
    for repeat in range(n_repeats):
        # get simulated causal genes
        sim_causal_gene_file = f"{exp_name}/sim_repeat_{repeat}/simulation_data/causal_genes_simulated.parquet"
        sim_causal_genes = pd.read_parquet(sim_causal_gene_file)

        # get number of variants with given annotation in simulated causal genes
        vars_causal_genes = vars_to_ensgid_mapping[
            vars_to_ensgid_mapping["gene"].isin(sim_causal_genes["gene"])
        ]
        anno_frequency_gene_variants = annotation_df.loc[
            vars_causal_genes["variant_id"]
        ].sum(axis=0)

        # get simulated causal variants
        sim_causal_variant_file = f"{exp_name}/sim_repeat_{repeat}/simulation_data/simulated_causal_variants.pkl"
        with open(sim_causal_variant_file, "rb") as f:
            sim_causal_variants = pickle.load(f)

        vars_this_repeat = annotation_df.copy().loc[sim_causal_variants]
        vars_this_repeat = vars_this_repeat.sum(axis=0).to_frame().reset_index()
        vars_this_repeat.columns = ["consequence", "n_causal_variants"]
        vars_this_repeat = vars_this_repeat.assign(
            total_n_variants=lambda df: df["consequence"].map(
                lambda consequence: anno_frequency_gene_variants.get(consequence, 0)
            )
        )

        # relative count = number of causal variants with a given annotation divided by the total number of variants
        # with that annotation in all simulated causal genes
        vars_this_repeat["relative_count"] = (
            vars_this_repeat["n_causal_variants"] / vars_this_repeat["total_n_variants"]
        )
        causal_vars_with_anno_list.append(vars_this_repeat)

        ## Number of causal variants per gene
        n_causal_vars_per_gene = (
            vars_causal_genes.query("variant_id in @sim_causal_variants")
            .drop(columns="gene_base")
            .groupby("gene")
            .size()
            .reset_index()
            .rename(columns={0: "n_causal_variants"})
        )
        causal_vars_per_gene_list.append(n_causal_vars_per_gene)

    causal_vars_with_anno = (
        pd.concat(causal_vars_with_anno_list)
        .groupby("consequence")
        .mean()
        .reset_index()
    )
    n_causal_vars_per_gene = (
        pd.concat(causal_vars_per_gene_list).groupby("gene").mean().reset_index()
    )

    return causal_vars_with_anno, n_causal_vars_per_gene


def _get_data_for_variant_annoation(config):
    vars_to_ensgid_mapping = pd.read_parquet(
        config["simulation"]["simulation_config"]["vars_to_ensgid_mapping_file"]
    )
    vars_to_ensgid_mapping["variant_id"] = vars_to_ensgid_mapping["variant_id"].astype(
        int
    )

    anno_columns = copy.deepcopy(
        config["simulation"]["data"]["dataset_config"]["annotations"]
    )
    cols_to_exclude = [
        "UKB_AF",
        "sift_score",
        "condel_score",
        "polyphen_score",
        "CADD_raw",
        "AbSplice_DNA",
        "PrimateAI_score",
    ]
    deepsea_cols = [f"DeepSEA_PC_{i}" for i in range(10)]
    cols_to_exclude = cols_to_exclude + deepsea_cols
    logger.info(f"Columns to remove for plotting: {cols_to_exclude}")
    for col in cols_to_exclude:
        try:
            anno_columns.remove(col)
        except:
            logger.info(f"Column {col} not in columns")
    anno_columns.append("CADD_PHRED")
    anno_columns.append("id")

    annotation_file = config["simulation"]["data"]["dataset_config"]["annotation_file"]
    annotation_df = pd.read_parquet(annotation_file, columns=anno_columns).set_index(
        "id"
    )
    vep_cols = [col for col in annotation_df if "Consequence_" in col]
    annotation_df["any_consquence"] = annotation_df[vep_cols].sum(axis=1)
    annotation_df = (
        annotation_df.assign(
            any_consquence=lambda df: df["any_consquence"].map(
                lambda x: 1 if x > 0 else 0
            )
        )
        .assign(
            CADD_PHRED_5=lambda df: df["CADD_PHRED"].map(lambda x: 1 if x > 5 else 0)
        )
        .drop(columns="CADD_PHRED")
    )

    anno_frequency = annotation_df.sum(axis=0)
    anno_frequency = anno_frequency.sort_values()

    return vars_to_ensgid_mapping, annotation_df, anno_frequency


def _plot_causal_variant_properties(causal_vars_with_anno, exp_name, limits):
    plots = []

    def x_labeller(labels):
        return [re.sub(r"^Consequence_|_variant", "", label) for label in labels]

    for y in ["n_causal_variants", "relative_count"]:
        p = (
            ggplot(
                causal_vars_with_anno, aes(y=y, x="consequence")
            )  # , color = 'rare_maf'))
            + geom_bar(stat="identity", position="dodge", alpha=0.8)
            + theme_bw()
            + scale_x_discrete(labels=x_labeller, limits=limits)
            + theme(axis_text_x=element_text(angle=45, vjust=1.0, hjust=1.0))
        )
        # if y == 'n_causal_variants':
        #     n_causal_variants = int(re.search(r'n_total_causal_variants:(\d+)', exp_name).groups()[0])
        #     p = p + geom_hline(yintercept = n_causal_variants, color = 'darkred',
        #                        linetype = 'dashed', size = 0.8)
        plots.append(p)

    return plots


def _plot_number_of_vars_per_gene(n_causal_vars_per_gene):
    plots = []
    p_bar = p = (
        ggplot(
            n_causal_vars_per_gene, aes(y="n_causal_variants", x="gene")
        )  # , color = 'rare_maf'))
        + geom_bar(stat="identity", position="dodge", alpha=0.8)
        + ggtitle("Number of causal \n variants per gene \n (Average across repeats)")
        + theme_bw()
        # + scale_x_discrete(labels = x_labeller, limits = limits)
        + theme(axis_text_x=element_text(angle=45, vjust=1.0, hjust=1.0))
    )
    plots.append(p_bar)

    p_box = (
        ggplot(
            n_causal_vars_per_gene, aes(y="n_causal_variants", x="0")
        )  # , color = 'rare_maf'))
        + geom_violin()
        + geom_boxplot(width=0.2)
        + ggtitle("Number of causal \n variants per gene \n (Average across repeats)")
        + theme_bw()
        + theme(
            axis_text_x=element_blank(),
            axis_title_x=element_blank(),
            figure_size=(4, 6),
        )
    )
    plots.append(p_box)

    return plots


@click.group()
def cli():
    pass


@cli.command()
@click.option("--n-repeats", type=int, default=1)
@click.argument("exp-name", type=click.Path(exists=True))
@click.argument("config-file", type=click.Path(exists=True))
# @click.argument('sim-causal-genes', type=click.Path(exists=True))
def evaluate(n_repeats, exp_name, config_file):
    with open(config_file) as f:
        config = yaml.safe_load(f)

    out_path = f"{exp_name}/eval"
    b_vtypes = config["baseline"].get("variant_types", ["plof"])
    b_ttypes = config["baseline"].get("test_types", ["burden"])
    baseline_test_combis = [combi for combi in itertools.product(b_vtypes, b_ttypes)]
    logger.info(baseline_test_combis)

    metrics, all_associations = evaluate_(exp_name, n_repeats, baseline_test_combis)

    logger.info("Saving results")
    with open(f"{out_path}/metrics.pkl", "wb") as f:
        pickle.dump(metrics, f)
    all_associations.to_parquet(f"{out_path}/all_associations.parquet")


@cli.command()
@click.option("--n-repeats", type=int, default=1)
@click.argument("exp-name", type=str)
@click.argument("config-file", type=click.Path(exists=True))
def analyse_variant_properties(config_file, exp_name, n_repeats):
    with open(config_file) as f:
        config = yaml.safe_load(f)

    out_path = Path(f"{exp_name}/eval")
    logger.info("Loading required data")
    (
        vars_to_ensgid_mapping,
        annotation_df,
        anno_frequency,
    ) = _get_data_for_variant_annoation(config)

    logger.info("Retrieving causal variant properties")
    causal_vars_with_anno, n_causal_vars_per_gene = _get_causal_variant_properties(
        exp_name, n_repeats, annotation_df, vars_to_ensgid_mapping
    )
    plots = {}
    plots_var_properties = _plot_causal_variant_properties(
        causal_vars_with_anno, exp_name, list(anno_frequency.index)
    )
    plots["causal_variant_properties"] = plots_var_properties

    plots_variants_per_gene = _plot_number_of_vars_per_gene(n_causal_vars_per_gene)
    plots["variants_per_gene"] = plots_variants_per_gene

    logger.info("saving plots")
    for name, p in plots.items():
        if len(p) > 0:
            save_as_pdf_pages(p, filename=str(out_path / f"{name}.pdf"))
        else:
            (out_path / f"{name}.pdf").touch()


if __name__ == "__main__":
    cli()
