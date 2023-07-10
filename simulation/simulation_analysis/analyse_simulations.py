from sim_evaulation_functions import collect_all_results, get_maf_filter_var_explained
import sys
import os
import logging
import math
import itertools
import pickle
import pandas as pd
import numpy as np
import yaml

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


####### PATHS to be defined by user! ######

config_file = "config_eval.yaml"
with open(config_file) as f:
    config_eval = yaml.safe_load(f)
sim_res_dir = config_eval["sim_res_dir"]

r_plot_data_path = config_eval["r_plot_data_path"]
default_config = config_eval["default_config"]
############################
DEEPRVAT_ANALYSIS_DIR = os.environ["DEEPRVAT_ANALYSIS_DIR"]
with open(f"{DEEPRVAT_ANALYSIS_DIR}/data/correlated_genes_to_remove.pickle", "rb") as f:
    genes_to_exclude_dict = pickle.load(f)

tests_to_keep = [
    "Burden pLOF",
    "Burden missense",
    "SKAT pLOF",
    "SKAT missense",
    "Burden/SKAT combined",
    "DeepRVAT",
]

default_config_dir = default_config.replace(";", "_")


quantiles = [
    0.5
]  # how many strongly correlated genes should be excluded from the analysis
# strongly correlated genes lead to inflation of the true-positives for the baseline/seed gene discovery methods

################### Evaluate MAF filter experiments (Fig. 2b) ##############################
##########################################################################################


logger.info("Retrieving results for MAF filter experiments")

metrics_assoc_maf_1, avg_assoc_maf_1 = collect_all_results(
    exp_name=config_eval["exp_name_assoc_maf_1"],
    quantiles=quantiles,
    local_res_dir=sim_res_dir,
    config_cols=[
        "vprop",
        "maxaf",
        "vgen",
        "vnoise",
        "vvarnoise",
        "vcont",
        "vbin",
        "vmaf",
        "expbin",
    ],
    genes_to_exclude=genes_to_exclude_dict,
)

metrics_assoc_maf_01, avg_assoc_maf_01 = collect_all_results(
    exp_name=config_eval["exp_name_assoc_maf_01"],
    quantiles=quantiles,
    local_res_dir=sim_res_dir,
    config_cols=[
        "vprop",
        "maxaf",
        "vgen",
        "vnoise",
        "vvarnoise",
        "vcont",
        "vbin",
        "vmaf",
        "expbin",
    ],
    genes_to_exclude=genes_to_exclude_dict,
)
metrics_assoc_maf_001, avg_assoc_maf_001 = collect_all_results(
    exp_name=config_eval["exp_name_assoc_maf_001"],
    quantiles=quantiles,
    local_res_dir=sim_res_dir,
    config_cols=[
        "vprop",
        "maxaf",
        "vgen",
        "vnoise",
        "vvarnoise",
        "vcont",
        "vbin",
        "vmaf",
        "expbin",
    ],
    genes_to_exclude=genes_to_exclude_dict,
)


maf_categories = ["1%", "0.1%", "0.01%"]
avg_metrics_maf_filter_combined = (
    pd.concat(
        [avg_assoc_maf_1, avg_assoc_maf_01, avg_assoc_maf_001], keys=maf_categories
    )
    .reset_index(level=0)
    .rename(columns={"level_0": "assoc_maf"})
)

assoc_maf_cat = pd.Categorical(
    avg_metrics_maf_filter_combined["assoc_maf"], categories=maf_categories
)
avg_metrics_maf_filter_combined = avg_metrics_maf_filter_combined.assign(
    assoc_maf=assoc_maf_cat
)


test_type_cat = pd.Categorical(
    avg_metrics_maf_filter_combined["test"], categories=tests_to_keep
)
avg_metrics_maf_filter_combined = avg_metrics_maf_filter_combined.assign(
    test=test_type_cat
)


avg_metrics_maf_filter_combined.to_parquet(
    f"{r_plot_data_path}/avg_metrics_maf_filter_combined.parquet"
)


# ##################### Square cooredinate for plotting ###########################################

maf_categories = [0.01, 0.001, 0.0001]


def set_outline_color(x, y):
    if x == y:
        return "matched"
    else:
        return "unmatched"


combis = [
    [0.01, 0.01],
    [0.01, 0.001],
    [0.01, 0.0001],
    [0.001, 0.01],
    [0.001, 0.001],
    [0.001, 0.0001],
    [0.0001, 0.01],
    [0.0001, 0.001],
    [0.0001, 0.0001],
]
outline = pd.DataFrame(columns=["max_rare_af", "assoc_maf"], data=combis)
outline["outline_color"] = outline.apply(
    lambda x: set_outline_color(x["max_rare_af"], x["assoc_maf"]), axis=1
)
outline["outline_size"] = outline.apply(
    lambda x: set_outline_color(x["max_rare_af"], x["assoc_maf"]), axis=1
)

square = pd.DataFrame(
    {
        "x": [-math.inf, math.inf, math.inf, -math.inf],
        "y": [-math.inf, -math.inf, math.inf, math.inf],
    }
)

squere_outline = outline.merge(square, how="cross")
assoc_maf_cat = pd.Categorical(squere_outline["assoc_maf"], categories=maf_categories)
squere_outline = squere_outline.assign(assoc_maf=assoc_maf_cat)
squere_outline.to_parquet(f"{r_plot_data_path}/square_outline_maf_plot.parquet")

##################### Variance explained ###########################################
var_weights_per_bin = get_maf_filter_var_explained(
    f"{sim_res_dir}/{config_eval['exp_name_assoc_maf_1']}"
)
# can be any of the maf filter experiments because the data is the same
var_weights_per_bin.to_parquet(f"{r_plot_data_path}/var_weights_per_bin.parquet")


# ################### Evaluate variation of prop causal variants (Supp. Fig. 2.3a) ##############################
# ##################################################################################################################
logger.info("Retrieving results for prop causal variant experiments")

(
    metrics_vary_prop_causal,
    avg_vary_prop_causal,
    assoc_vary_prop_causal,
) = collect_all_results(
    exp_name=config_eval["exp_name_vary_prop_causal"],
    quantiles=quantiles,
    local_res_dir=sim_res_dir,
    return_associations=True,
    genes_to_exclude=genes_to_exclude_dict,
)
test_type_cat = pd.Categorical(avg_vary_prop_causal["test"], categories=tests_to_keep)
avg_vary_prop_causal = avg_vary_prop_causal.assign(test=test_type_cat)

avg_vary_prop_causal.to_parquet(f"{r_plot_data_path}/avg_vary_prop_causal.parquet")


# ################### Replication plot (Supp. Fig. 2.3b) ##############################
# #######################################################################################


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
dfs_for_replication_plot = []
default_config_assoc_result = assoc_vary_prop_causal.query("config == @default_config")
for test_name, tests in BASELINE_METHODS_DICT.items():
    this_baseline_df = (
        assoc_vary_prop_causal.copy()
        .query('combination == "Baseline" & method in @tests')
        .assign(Method=test_name)
    )
    dfs_for_replication_plot.append(this_baseline_df)

deeprvat_df = (
    assoc_vary_prop_causal.copy()
    .query('combination == "DeepRVAT"')
    .assign(Method="DeepRVAT")
)
dfs_for_replication_plot.append(deeprvat_df)
replication_results = pd.concat(dfs_for_replication_plot)
replication_results["replicated"] = replication_results["causal_gene"]


n_genes_per_pheno = 200
top_list = []
for name, group in replication_results.groupby(["repeat", "Method"]):
    this_df = group.sort_values("pval_corrected")
    this_df = this_df.drop_duplicates(subset=["gene"])
    this_df = this_df.head(n_genes_per_pheno)
    this_df["Gene rank"] = np.arange(1, n_genes_per_pheno + 1)
    this_df["Replicated genes"] = np.cumsum(this_df["replicated"])
    if name[1].startswith("deeprvat"):
        this_df["Method"] = "deeprvat"
    top_list.append(this_df)

top = pd.concat(top_list)

top["Significant"] = pd.Categorical(
    top["significant"], categories=[True, False], ordered=True
)
top["Method"] = pd.Categorical(
    top["Method"],
    categories=[
        "Burden pLOF",
        "Burden missense",
        "SKAT pLOF",
        "SKAT missense",
        "Burden/SKAT combined",
        "DeepRVAT",
    ],
    ordered=True,
)

top_with_uncertainty = (
    top[["Gene rank", "Replicated genes", "Method", "repeat"]]
    .groupby(["Method", "Gene rank"])
    .agg({"Replicated genes": ["mean", "std"]})
    .reset_index()[["Method", "Gene rank", "Replicated genes"]]
    .set_index(["Method", "Gene rank"])["Replicated genes"]
    .reset_index()
)
top_with_uncertainty
top_with_uncertainty.to_parquet(f"{r_plot_data_path}/top_with_uncertainty.parquet")


################### FDR calibration (Supp. Fig. 2.2) ############################################################
##################################################################################################################

logger.info("Retrieving data for FDDR calibration plots (runs long)")
intervals = list(np.linspace(0.01, 0.1, 20))
for i in [0.05, 0.1, 0.2]:
    intervals.append(i)

alphas_to_test = intervals
fdr_metrics_list = []
fdr_avg_metrics_list = []

fdr_data_dir = "./fdr_data/"
if not os.path.exists(fdr_data_dir):
    os.makedirs(fdr_data_dir)
for alpha in alphas_to_test:
    print(f"\n {alpha} \n")
    fdr_metrics_vary_prop_causal, fdr_avg_vary_prop_causal = collect_all_results(
        exp_name=config_eval["exp_name_vary_prop_causal"],
        quantiles=[0.5],
        local_res_dir=sim_res_dir,
        alpha=alpha,
        genes_to_exclude=genes_to_exclude_dict,
    )

    fdr_metrics_list.append(fdr_metrics_vary_prop_causal)
    fdr_avg_metrics_list.append(fdr_avg_vary_prop_causal)

    fdr_res = {
        "metrics": fdr_metrics_vary_prop_causal,
        "avg_metrics": fdr_avg_vary_prop_causal,
    }
    with open(f"{fdr_data_dir}/fdr_{alpha}_plot_data.pickle", "wb") as f:
        pickle.dump(fdr_res, f)

fdr_avg_metrics = pd.concat(
    [
        fdr_avg_metrics_list[i].assign(expected_FDR=alphas_to_test[i])
        for i in range(len(fdr_avg_metrics_list))
    ]
).rename(columns={"expected_FDR": "expected FDR"})

fdr_avg_metrics.to_parquet(f"{r_plot_data_path}/fdr_avg_metrics.parquet")
