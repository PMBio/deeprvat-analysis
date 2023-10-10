import copy
import logging
import pickle
import sys
from pathlib import Path
from pprint import pprint

import click
import pandas as pd
import numpy as np
import plotnine as p9
import yaml
from deeprvat.utils import pval_correction
from tqdm import tqdm

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)



PHENOTYPES = [
    "Apolipoprotein_A",
    "Apolipoprotein_B",
    "Calcium",
    "Cholesterol",
    "HDL_cholesterol",
    "IGF_1",
    "LDL_direct",
    "SHBG",
    "Total_bilirubin",
    "Triglycerides",
    "Urate",
    "Mean_corpuscular_volume",
    "Platelet_count",
    "Mean_platelet_thrombocyte_volume",
    "Platelet_crit",
    "Standing_height",
    "Mean_reticulocyte_volume",
    "Platelet_distribution_width",
    "Lymphocyte_percentage",
    "Neutrophill_count",
    "Red_blood_cell_erythrocyte_count",
]




METHODS = list(
    reversed(
        [f"{m} {t}" for m in ["Burden", "SKAT"] for t in ["pLOF", "missense"]]
        + [  # , '"likely deleterious"']] +
            "Burden pLOF/missense",
            "SKAT pLOF/missense",
            "Burden/SKAT combined",
            "DeepRVAT",
        ]
    )
)


def read_comparison_results(comparison_dir: str, gene_df: pd.DataFrame, pheno: str):
    pheno = 'IGF_1'
    phenocode_df = pd.read_parquet("phenocodes.parquet", engine="pyarrow")
    pheno = 'IGF-1' if pheno == 'IGF_1' else pheno
    phenocode_dict = {
        pheno: int(code.split("-")[0]) if code.split("-")[0].isdigit() else code
        for pheno, code in zip(phenocode_df["phenotype"], phenocode_df["phenocode"])
    }

    GENEBASS_NAME_DICT = copy.deepcopy(phenocode_dict)
    # GENEBASS_NAME_DICT.update({f"{k}_standardized": v for k, v in phenocode_dict.items()})

    UKB_500_NAME_DICT = {
        pheno: trait
        for pheno, trait in zip(
            phenocode_df["phenotype"], phenocode_df["backman_trait"]
        )
    }
    UKB_500_NAME_DICT.update(
        {
            pheno: trait
            for pheno, trait in zip(
                phenocode_df["phenotype"], phenocode_df["backman_trait"]
            )
        }
    )

    gene_df = gene_df.copy()
    gene_df["ensgid"] = gene_df["gene"].str.split(".", expand=True)[0]

    genebass = pd.read_parquet(
        Path(comparison_dir) / "genebass_pvals_500k_selected_traits.parquet",
        engine="pyarrow",
    )
    genebass_phenocode = str(GENEBASS_NAME_DICT[pheno])
    genebass = genebass.query(
        "phenocode == @genebass_phenocode"
        " and ((significant_burden and keep_gene_burden) or (significant_skato and keep_gene_skato))"
    )
    # significant_mask = (genebass[[
    #     c for c in genebass.columns if c.startswith('significant_')
    # ]].sum(axis=1) > 0)
    # genebass_ensgids = genebass.loc[significant_mask, 'gene_id']
    genebass_ensgids = genebass["gene_id"]
    genebass_ids = gene_df.query("ensgid in @genebass_ensgids")["id"]

    backman = pd.read_excel(
        Path(comparison_dir) / "41586_2021_4103_MOESM5_ESM.xlsx",
        sheet_name="SD2",
        engine="openpyxl",
    )
    backman = backman[backman["Marker type"] == "Burden"]
    backman_pheno = UKB_500_NAME_DICT[pheno]
    markers = ["M1.1", "M3.1"]
    backman = backman.query("Trait == @backman_pheno and Marker in @markers")
    backman_genenames = backman["Gene"]
    backman_ids = gene_df.query("gene_name in @backman_genenames")["id"]

    return {
        "genebass": genebass_ids,
        "UKB500k": backman_ids,
    }

def prep_for_rep_plot(
    plotting_results,
    comparison_results,
    n_genes=600,
    n_genes_per_pheno=100,
    title=None,
    method_mapping=None,
    phenotype_mapping=None,
    colors=None,
    skip_pheno_plots=True,
):
    replication_results = plotting_results.copy()
    replication_results = replication_results.rename(
        columns={"experiment_group": "Method", "significant": "Significant"}
    )

    if method_mapping is not None:
        replication_results["Method"] = replication_results["Method"].apply(
            lambda x: method_mapping[x] if x in method_mapping else x
        )

    replication_results = replication_results.query("Method in @METHODS")
    replication_results["Method"] = pd.Categorical(
        replication_results["Method"], categories=METHODS, ordered=True
    )

    # replication_results = results.query(f'((experiment_group == @selected_experiment and experiment == "DeepRVAT ({n_repeats} repeats)") or '
    #                                    '(experiment_group in @baseline_groups)) '
    #                                    ' and correction_method == "FDR"'
    #                                    ' and phenotype in @pheno_to_use'
    #                                   )

    replication_results["replicated"] = [
        gene in comparison_results[pheno]
        for pheno, gene in zip(
            replication_results["phenotype"], replication_results["gene"]
        )
    ]

    if phenotype_mapping is not None:
        replication_results["phenotype"] = replication_results["phenotype"].apply(
            lambda x: phenotype_mapping[x]
            if x in phenotype_mapping
            else " ".join(x.replace("_standardized", "").split("_"))
        )

    top_list = []
    for name, group in replication_results.groupby("Method"):
        this_df = group.sort_values("pval_corrected")
        this_df = this_df.drop_duplicates(subset=["phenotype", "gene"])
        this_df = this_df.head(n_genes)
        this_df["Gene rank"] = np.arange(1, n_genes + 1)
        this_df["Replicated genes"] = np.cumsum(this_df["replicated"])
        if name.startswith("deeprvat"):
            this_df["Method"] = "deeprvat"
        top_list.append(this_df)

    top = pd.concat(top_list)
    top["Significant"] = pd.Categorical(
        top["Significant"], categories=[True, False], ordered=True
    )
    top = top[
        [
            "Gene rank",
            "Replicated genes",
            "Method",
            "Significant",
            "phenotype",
            "replicated",
            "Discovery type",
            "gene",
        ]
    ].rename(columns={"phenotype": "Trait"})
    return top


@click.command()
@click.option("--out-dir", type=click.Path(exists=True), default=".")
@click.option("--recompute-comparison-results", is_flag=True)
@click.option("--analyze-all-repeats", is_flag=True)
@click.argument("experiment-dir", type=click.Path(exists=True))
def cli(
        out_dir: str, 
        experiment_dir: str, 
        recompute_comparison_results: bool,
        analyze_all_repeats: bool,
):
    logger.info(f'analysing phenotypes {PHENOTYPES}')

    if recompute_comparison_results:
        comparison_results = {
            pheno: read_comparison_results(
                ".", pd.read_parquet("protein_coding_genes.parquet"), pheno
            )
            for pheno in tqdm(PHENOTYPES)
        }
        comparison_results = {
            " ".join(pheno.split("_")): set(pd.concat(r.values()))
            for pheno, r in comparison_results.items()
        }
        with open("comparison_results.pkl", "wb") as f:
            pickle.dump(comparison_results, f)
         df = pd.DataFrame([(key, value) for key, values in comparison_results.items() for value in values], columns=['phenotype', 'gene'])
        df['phenotype'] = [name.replace(' ', '_') for name in df['phenotype']]
        df.to_parquet('comparison_results.parquet') #deeprvat-analysis/data/comparison_results.parquet" used by monti/staar replication scripts
    else:
        with open("comparison_results.pkl", "rb") as f:
            comparison_results = pickle.load(f)

    results = pd.concat(
        [
            pd.read_parquet(
                Path(experiment_dir) / p / "deeprvat/eval/all_results.parquet",
                engine="pyarrow",
            )
            for p in PHENOTYPES
        ]
    )
    phenotypes_to_remove = set(results['phenotype'].unique()) - set(comparison_results.keys())
    logger.info(f'excluding pheotypes {phenotypes_to_remove} because they are not in comparison studies')
    results = results[~results['phenotype'].isin(phenotypes_to_remove)]

    with open(Path(experiment_dir) / "config.yaml") as f:
        config = yaml.safe_load(f)
    n_repeats = config["n_repeats"]

    if analyze_all_repeats:
        repeats_to_use = range(1, n_repeats + 1)
    else:
        repeats_to_use = [n_repeats]

    all_repeats_list = []
    for repeats in repeats_to_use:
        logger.info(f"Analyzing replication with {repeats} DeepRVAT repeats")
        rep_list = []
        replication_data_all = prep_for_rep_plot(
            results.query("repeats == @repeats"),
            comparison_results,
            n_genes=1000,
        ).assign(pheno_grouping="all_phenotypes")
        rep_list.append(replication_data_all)

        for pheno in results["phenotype"].unique():
            replication_data_pheno = prep_for_rep_plot(
                results.query("phenotype == @pheno"),
                comparison_results,
                n_genes=1000,
            ).assign(pheno_grouping="single_pheno")
            rep_list.append(replication_data_pheno)

        this_replication_data = pd.concat(rep_list).query("Method in @METHODS")
        this_replication_data["repeats"] = repeats
        all_repeats_list.append(this_replication_data)

    logger.info("Writing replication")
    replication_data = pd.concat(all_repeats_list)
    replication_data.to_parquet(Path(out_dir) / "replication.parquet", engine="pyarrow")


if __name__ == "__main__":
    cli()
