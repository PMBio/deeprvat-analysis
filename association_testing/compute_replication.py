import copy
import pickle
from pathlib import Path
from pprint import pprint

import click
import pandas as pd
import numpy as np
import plotnine as p9
import yaml
from deeprvat.utils import pval_correction
from tqdm import tqdm

phenocode_df = pd.read_parquet("phenocodes.parquet", engine="pyarrow")
# phenocode_dict = {
#     pheno: int(code.split("-")[0])
#     for pheno, code in zip(phenocode_df["phenotype"], phenocode_df["phenocode"])
# }
phenocode_dict = {}
for pheno, code in zip(phenocode_df['phenotype'], phenocode_df['phenocode']):
    try:
        code = int(code.split('-')[0])
    except:
        code = code
    phenocode_dict[pheno] = code


GENEBASS_NAME_DICT = copy.deepcopy(phenocode_dict)
# GENEBASS_NAME_DICT.update({f"{k}_standardized": v for k, v in phenocode_dict.items()})

UKB_500_NAME_DICT = {
    pheno: trait
    for pheno, trait in zip(phenocode_df["phenotype"], phenocode_df["backman_trait"])
}
UKB_500_NAME_DICT.update(
    {
        pheno: trait
        for pheno, trait in zip(
            phenocode_df["phenotype"], phenocode_df["backman_trait"]
        )
    }
)

OLD_PHENOTYPES = [
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
TRAINING_PHENOTYPES = copy.deepcopy(OLD_PHENOTYPES)

# PHENOTYPES = [
#     "Apolipoprotein_A",
#     "Apolipoprotein_B",
#     "Calcium",
#     "Cholesterol",
#     "HDL_cholesterol",
#     "IGF_1",
#     "LDL_direct",
#     "SHBG",
#     "Total_bilirubin",
#     "Triglycerides",
#     "Urate",
#     "Mean_corpuscular_volume",
#     "Platelet_count",
#     "Mean_platelet_thrombocyte_volume",
#     "Platelet_crit",
#     "Standing_height",
#     "Mean_reticulocyte_volume",
#     "Platelet_distribution_width",
#     "Lymphocyte_percentage",
#     "Neutrophill_count",
#     "Red_blood_cell_erythrocyte_count",
#     "Body_mass_index_BMI",
#     "Glucose",
#     "Vitamin_D",
#     "Albumin",
#     "Total_protein",
#     "Cystatin_C",
#     "Gamma_glutamyltransferase",
#     "Alkaline_phosphatase",
#     "Creatinine"
# ]

NEW_PHENOTYPES = [
    "Body_mass_index_BMI",
    "Glucose",
    "Vitamin_D",
    "Albumin",
    "Total_protein",
    "Cystatin_C",
    "Gamma_glutamyltransferase",
    "Alkaline_phosphatase",
    "Creatinine", 
    "Whole_body_fat_free_mass",
    "Forced_expiratory_volume_in_1_second_FEV1",
    "Glycated_haemoglobin_HbA1c",
    "QTC_interval",
    'WHR',
    'WHR_Body_mass_index_BMI_corrected' #not using these phenos because they are not in genebas/backman
]

# PHENOTYPES = NEW_PHENOTYPES
PHENOTYPES = [*NEW_PHENOTYPES, *OLD_PHENOTYPES]
# PHENOTYPES = OLD_PHENOTYPES
print(f'analysing phenotypes {PHENOTYPES}')

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
@click.argument("experiment-dir", type=click.Path(exists=True))
def cli(out_dir: str, experiment_dir: str, recompute_comparison_results: bool):

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
    else:
        with open("comparison_results.pkl", "rb") as f:
            comparison_results = pickle.load(f)

    result_list = []
    for p in PHENOTYPES:
        this_res = pd.read_parquet(
               
                Path(experiment_dir) / p / "deeprvat/eval/all_results.parquet",
                engine="pyarrow",
            )
        if p in NEW_PHENOTYPES:
            print(f'using DeepRVAT wo baseline counts for DeepRVAT for pheno {p}')
            this_res = this_res.query('Method != "DeepRVAT"')
            this_res.loc[this_res['Method'] == 'DeepRVAT wo baseline', 'Method'] = 'DeepRVAT'
        result_list.append(this_res)
    results = pd.concat(result_list)
    phenotypes_to_remove = set(results['phenotype'].unique()) - set(comparison_results.keys())
    print(f'excluding pheotypes {phenotypes_to_remove} because they are not in comparison studies')
    results = results[~results['phenotype'].isin(phenotypes_to_remove)]
    # results = pd.concat(
    #     [
    #         pd.read_parquet(
    #             Path(experiment_dir) / p / "deeprvat/eval/all_results.parquet",
    #             engine="pyarrow",
    #         )
    #         for p in PHENOTYPES
    #     ]
    # )
    rep_list = []
    replication_data_all = prep_for_rep_plot(
        results,
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

    replication_data = pd.concat(rep_list).query("Method in @METHODS")

    print("Writing replication")
    replication_data.to_parquet(Path(out_dir) / "replication.parquet", engine="pyarrow")


if __name__ == "__main__":
    cli()
