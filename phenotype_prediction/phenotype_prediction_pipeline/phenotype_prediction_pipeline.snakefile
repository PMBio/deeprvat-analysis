from snakemake.utils import Paramspace
from snakemake.utils import Paramspace
from snakemake.utils import min_version
import copy
import os

min_version("6.0")


configfile: "config_eval.yaml"


from pathlib import Path
import pandas as pd

conda_check = 'conda info | grep "active environment"'

phenotypes = config.get("phenotypes")


DEEPRVAT_ANALYSIS_DIR = os.environ["DEEPRVAT_ANALYSIS_DIR"]
code_dir = (
    f"{DEEPRVAT_ANALYSIS_DIR}/phenotype_prediction/phenotype_prediction_pipeline/"
)

top_bottom_quantiles = ["topq", "bottomq"]

regression_config = config["regression"]
regression_config["r_config"]["code_dir"] = code_dir
phenotype_suffixes = regression_config["pheno_suffixes"]
top_quantiles = regression_config["top_quantiles"]

btypes = copy.deepcopy(regression_config["r_config"]["btypes"])
print(btypes)
# btypes.append("deeprvat")
# btypes.remove('sift')
burden_btypes = [btype for btype in btypes if "_" not in btype]

py = f"python {code_dir}/"

phenotypes = config['phenotypes']



# To run the phenotype prediction pipeline only (i.e., get predictions for each test sample and compute metrics)
# only run the pipeline using `rule all_data`.
# To also get all required data for the paper plots, run the complete snakemake workflow.


deeprvat_run_dir = config.get('deeprvat_run_dir')
#TODO update this
burden_dir_mapper = {

    btype: f"alternative_burdens/{btype}" if btype != "deeprvat" else "deeprvat_burdens"
    for btype in burden_btypes
}

burden_phenotype = phenotypes[0]



rule all_phenotype_prediction:
    input:
        expand(
            "models/logistic_models/plotting_data/combined_metrics_{phenotype_suffix}.Rds",
            phenotype_suffix=phenotype_suffixes,
        ),
        expand(
            "models/linear_models/plotting_data/all_recomputed_metrics_test_{phenotype_suffix}.Rds",
            phenotype_suffix=phenotype_suffixes,

        ),
        expand(
            "models/linear_models/plotting_data/replication_in_extreme/ranked_df_combined_{outlier_mode}_{extreme_quantile}_{phenotype_suffix}.Rds",
            phenotype_suffix=phenotype_suffixes,
            outlier_mode=["both", "bottomq", "topq"],
            extreme_quantile=[0.01, 0.001, 0.005],
        ),
        expand(
            "models/linear_models/plotting_data/enrichment_prs_vs_zscore_{quantile}_{phenotype_suffix}.Rds",
            phenotype_suffix=phenotype_suffixes,
            quantile=[0.99, 0.999],

        ),



rule replication_in_extremes_plot:
    conda:
        "r-env"
    input:
        "models/linear_models/plotting_data/plot_df_list_sub_{phenotype_suffix}.finished",
    output:
        "models/linear_models/plotting_data/replication_in_extreme/ranked_df_combined_{outlier_mode}_{extreme_quantile}_{phenotype_suffix}.Rds",
    params:
        phenotype_suffix="{phenotype_suffix}",
        linear_model_res_path="models/linear_models",
        phenotypes=phenotypes,
        extreme_quantile="{extreme_quantile}",
        outlier_mode="{outlier_mode}",
    resources:
        mem_mb=40960,
        load=16000,
    script:
        f"{code_dir}/prep_paper_figure_data/replication_in_extremes_plot.R"


rule plot_enrichment_prs_vs_zscore:
    conda:
        "r-env"
    input:
        "models/linear_models/plotting_data/plot_df_list_sub_{phenotype_suffix}.finished",
    output:
        "models/linear_models/plotting_data/enrichment_prs_vs_zscore_{quantile}_{phenotype_suffix}.Rds",
    params:
        phenotype_suffix="{phenotype_suffix}",
        linear_model_res_path="models/linear_models",
        phenotypes=phenotypes,
        quantile="{quantile}",
    resources:
        mem_mb=40960,
        load=16000,
    script:
        f"{code_dir}/prep_paper_figure_data/plot_enrichment_prs_vs_zscore.R"


rule split_plot_df_list:
    conda:
        "r-env"
    input:
        "models/linear_models/plotting_data/plot_df_list_{phenotype_suffix}.Rds",
    output:
        "models/linear_models/plotting_data/plot_df_list_sub_{phenotype_suffix}.finished",
    params:
        phenotype_suffix="{phenotype_suffix}",
        linear_model_res_path="models/linear_models",
    resources:
        mem_mb=40960,
        load=16000,
    script:
        f"{code_dir}/prep_paper_figure_data/split_plot_df_list.R"


rule make_plot_df_list_data_linear_model:
    conda:
        "r-env"
    input:
        expand(
            "models/linear_models/{phenotype}_{{phenotype_suffix}}.Rds",
            phenotype=phenotypes,
        ),
    output:
        "models/linear_models/plotting_data/plot_df_list_{phenotype_suffix}.Rds",
    params:
        input_res_path="models/linear_models",
        out_dir="models/linear_models/plotting_data/",
        phenotype_suffix="{phenotype_suffix}",
        phenotypes=phenotypes,
    resources:
        mem_mb=40960,
        load=16000,
    script:
        f"{code_dir}/prep_paper_figure_data/make_plot_df_list_data.R"


rule all_data:
    input:
        expand(
            "models/linear_models/plotting_data/all_recomputed_metrics_test_{phenotype_suffix}.Rds",
            phenotype_suffix=phenotype_suffixes,
        ),
        expand(
            "models/logistic_models/plotting_data/combined_metrics_{phenotype_suffix}.Rds",
            phenotype_suffix=phenotype_suffixes,
        ),


rule prep_metrics_linear_model:
    conda:
        "r-env"
    input:
        expand(
            "models/linear_models/{phenotype}_{{phenotype_suffix}}.Rds",
            phenotype=phenotypes,
        ),
    output:
        "models/linear_models/plotting_data/all_recomputed_metrics_test_{phenotype_suffix}.Rds",
    params:
        input_res_path="models/linear_models/",
        out_dir="models/linear_models/plotting_data/",
        phenotype_suffix="{phenotype_suffix}",
        phenotypes=phenotypes,
    resources:
        mem_mb=20480,
        load=16000,
    script:
        f"{code_dir}/prep_linear_model_plotting_metrics.R"


rule prep_metrics_logistic_model:
    conda:
        "r-env"
    input:
        expand("models/logistic_models/{phenotype}_{{phenotype_suffix}}_{top_bottom_q}-{top_q}.Rds",
            phenotype=phenotypes,
            top_q=top_quantiles,
            top_bottom_q=top_bottom_quantiles,
        ),
    output:
        metrics="models/logistic_models/plotting_data/combined_metrics_{phenotype_suffix}.Rds",
    params:
        input_res_path="models/logistic_models/",
        out_dir="models/logistic_models/plotting_data/",
        phenotype_suffix="{phenotype_suffix}",
        phenotypes=phenotypes,
        top_quantiles=top_quantiles,
        top_bottom_q_vals=top_bottom_quantiles,
        code_dir=code_dir,
    resources:
        mem_mb=lambda wildcards, attempt: 60480 + (attempt - 1) * 10098,
        load=16000,
    script:
        f"{code_dir}/prep_logistic_model_plotting_metrics.R"



rule all_regression:
    input:
        expand(
            "models/linear_models/{phenotype}_{phenotype_suffix}.Rds",
            phenotype=phenotypes,
            phenotype_suffix=phenotype_suffixes,
        ),
        expand(
            "models/logistic_models/{phenotype}_{phenotype_suffix}_{top_bottom_q}-{top_q}.Rds",
            phenotype=phenotypes,
            phenotype_suffix=phenotype_suffixes,
            top_q=top_quantiles,
            top_bottom_q=top_bottom_quantiles,
        ),


rule fit_logistic_regression_model:
    conda:
        "r-env"
    input:
        "r_data/{phenotype}/data.finished",
    output:
        "models/logistic_models/{phenotype}_{phenotype_suffix}_{top_bottom_q}-{top_q}.Rds",
    resources:
        mem_mb=lambda wildcards, attempt: 20480 + (attempt - 1) * 4098,
        load=16000,
    params:
        model_out_file="models/logistic_models/{phenotype}_{phenotype_suffix}_{top_bottom_q}-{top_q}.Rds",
        phenotype="{phenotype}",
        phenotype_suffix="{phenotype_suffix}",
        top_q="{top_q}",
        regression_data_input_path="r_data/{phenotype}",
        config=regression_config["r_config"],
        top_bottom_q="{top_bottom_q}",
    script:
        f"{code_dir}/run_logistic_regression.R"


rule fit_linear_regression_model:
    conda:
        "r-env"
    input:
        "r_data/{phenotype}/data.finished",
    output:
        "models/linear_models/{phenotype}_{phenotype_suffix}.Rds",
    resources:
        mem_mb=lambda wildcards, attempt: 20480 + (attempt - 1) * 4098,
        load=16000,
    params:
        model_out_file="models/linear_models/{phenotype}_{phenotype_suffix}.Rds",
        phenotype="{phenotype}",
        phenotype_suffix="{phenotype_suffix}",
        regression_data_input_path="r_data/{phenotype}",
        config=regression_config["r_config"],
    script:
        f"{code_dir}/run_linear_regression.R"


rule all_r_input_data:
    input:
        expand(
            "r_data/{phenotype}/data.finished",
            phenotype=phenotypes,
        ),


rule prep_data_for_r:
    input:
        burdens=expand(
            "burdens/{btype}_burdens.pickle",
            btype=burden_btypes,
        ),
        config="config_eval.yaml",
        samples_test='train_test_samples/test_samples.pkl',
        samples_train='train_test_samples/train_samples.pkl',
        discoveries="discoveries/all_significant.parquet",
    output:
        "r_data/{phenotype}/data.finished",
    params:
        # dataset_dir="cv_split{cvsplit}/deeprvat/{phenotype}/deeprvat/",
        burdens_dir="burdens",
        out_dir="r_data/{phenotype}",
        burden_types=[f"--burden-type {b} " for b in burden_btypes],
    resources:
        mem_mb=64000,
        load=16000,
    shell:
        " && ".join(
        [
            conda_check,
            (
        py + "prep_data_for_regression_in_r.py prep-r-regression-data"
                    " --phenotype {wildcards.phenotype}"
                    " {params.burden_types}"
                    " --discoveries-file {input.discoveries}"
                    # " --dataset-dir {params.dataset_dir} "
                    " --burdens-dir {params.burdens_dir} "
                    " --out-dir {params.out_dir} "
                    " {input.samples_train}"
                    " {input.samples_test}"
                    " {input.config}"
                ),
                "touch {output}",
            ]
        )


rule all_burdens:
    input:
        expand(
            "burdens/{btype}_burdens.pickle",
            btype=burden_btypes,
        ),



# exctract burdens/ DeepRVAT gene impairment scores for significant genes
rule extract_burdens:
    input:
        discoveries=f"discoveries/all_significant.parquet",
        # burdens=lambda wildcards: f"cv_split{{cvsplit}}/{burden_dir_mapper[wildcards.btype]}/{burden_phenotype}/{wildcards.btype}/burdens{split_mapper[wildcards.split]}/burdens_zarr.created",
    output:
        "burdens/{btype}_burdens.pickle",
    resources:
        mem_mb=16000,
        load=16000,
    params:
        burden_input_dir=lambda wildcards: f"{burden_dir_mapper[wildcards.btype]}",
    shell:
        " && ".join(
        [
            conda_check,
            (
        py + "extract_burdens_for_prs.py extract-burdens "
                    " --burden-input-dir {params.burden_input_dir} "
                    " --burden-type {wildcards.btype}"
                    " --burdens-out {output}"
                    " --discoveries-file {input.discoveries}"
                ),
            ]
        )


rule all_discoveries:
    input:
        "discoveries/all_significant.parquet",


rule combine_discoveries:
    input:
        discoveries=expand(
            f"{deeprvat_run_dir}/{{phenotype}}/deeprvat/eval/significant.parquet",
            phenotype=phenotypes,
        ),
    output:
        "discoveries/all_significant.parquet",
    resources:
        mem_mb=8000,
        load=16000,
    params:
        sig_files=lambda wildcards, input: "".join(
            [f"--sig-file {d} " for d in input.discoveries]
        ),
    shell:
        " && ".join(
        [
            conda_check,
            (
        py + "extract_burdens_for_prs.py combine-significant "
                    " {params.sig_files} "
                    " --min-discoveries 2 "
                    " {output}"
                ),
            ]
        )


