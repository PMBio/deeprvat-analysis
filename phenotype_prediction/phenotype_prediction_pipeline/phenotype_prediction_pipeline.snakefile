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


cv_splits = config.get("cv_splits", 5)


phenotypes = config.get("phenotypes")
training_phenotypes = config.get("phenotypes")
phenotypes = ['HDL_cholesterol', 'LDL_direct']
training_phenotypes = phenotypes


DEEPRVAT_ANALYSIS_DIR = os.environ["DEEPRVAT_ANALYSIS_DIR"]
code_dir = (
    f"{DEEPRVAT_ANALYSIS_DIR}/phenotype_prediction/phenotype_prediction_pipeline/"
)

top_bottom_quantiles = ["topq", "bottomq"]

regression_config = config["regression"]
regression_config["r_config"]['code_dir'] = code_dir
phenotype_suffixes = regression_config["pheno_suffixes"]
top_quantiles = regression_config["top_quantiles"]

fdrs = regression_config.get("fdrs", [0.05])
btypes = copy.deepcopy(regression_config["r_config"]["vtypes"])
btypes.append("deeprvat")
burden_btypes = [btype for btype in btypes if "_" not in btype]

py = f"python {code_dir}/"


# To run the phenotype prediction pipeline only (i.e., get predictions for each test sample and compute metrics)
# only run the pipeline using `rule all_data`.
# To also get all required data for the paper plots, run the complete snakemake workflow.


rule all_phenotype_prediction:
    input:
        expand(
            "phenotype_prediction/models/linear_models/plotting_data/replication_in_extreme/ranked_df_combined_{outlier_mode}_{extreme_quantile}_{phenotype_suffix}_{fdr}.Rds",
            phenotype_suffix=phenotype_suffixes,
            outlier_mode=["both", "bottomq", "topq"],
            extreme_quantile=[0.01, 0.001, 0.005],
            fdr=fdrs,
        ),
        expand(
            "phenotype_prediction/models/linear_models/plotting_data/enrichment_prs_vs_zscore_{quantile}_{phenotype_suffix}_{fdr}.Rds",
            phenotype_suffix=phenotype_suffixes,
            quantile=[0.99, 0.999],
            fdr=fdrs,
        ),
        expand(
            "phenotype_prediction/models/linear_models/plotting_data/all_recomputed_metrics_test_{phenotype_suffix}_{fdr}.Rds",
            phenotype_suffix=phenotype_suffixes,
            fdr=fdrs,
        ),
        expand(
            "phenotype_prediction/models/logistic_models/plotting_data/combined_metrics_{phenotype_suffix}_{fdr}.Rds",
            phenotype_suffix=phenotype_suffixes,
            fdr=fdrs,
        ),


rule replication_in_extremes_plot:
    conda:
        "r-env"
    input:
        "phenotype_prediction/models/linear_models/plotting_data/plot_df_list_sub_{phenotype_suffix}_{fdr}.finished",
    output:
        "phenotype_prediction/models/linear_models/plotting_data/replication_in_extreme/ranked_df_combined_{outlier_mode}_{extreme_quantile}_{phenotype_suffix}_{fdr}.Rds",
    params:
        fdr="{fdr}",
        phenotype_suffix="{phenotype_suffix}",
        linear_model_res_path="phenotype_prediction/models/linear_models",
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
        "phenotype_prediction/models/linear_models/plotting_data/plot_df_list_sub_{phenotype_suffix}_{fdr}.finished",
    output:
        "phenotype_prediction/models/linear_models/plotting_data/enrichment_prs_vs_zscore_{quantile}_{phenotype_suffix}_{fdr}.Rds",
    params:
        fdr="{fdr}",
        phenotype_suffix="{phenotype_suffix}",
        linear_model_res_path="phenotype_prediction/models/linear_models/",
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
        "phenotype_prediction/models/linear_models/plotting_data/plot_df_list_{phenotype_suffix}_{fdr}.Rds",
    output:
        "phenotype_prediction/models/linear_models/plotting_data/plot_df_list_sub_{phenotype_suffix}_{fdr}.finished",
    params:
        fdr="{fdr}",
        phenotype_suffix="{phenotype_suffix}",
        linear_model_res_path="phenotype_prediction/models/linear_models/",
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
            "phenotype_prediction/models/linear_models/{phenotype}/{{phenotype_suffix}}_fdr-{{fdr}}.Rds",
            phenotype=phenotypes,
        ),
    output:
        "phenotype_prediction/models/linear_models/plotting_data/plot_df_list_{phenotype_suffix}_{fdr}.Rds",
    params:
        input_res_path="phenotype_prediction/models/linear_models/",
        out_dir="phenotype_prediction/models/linear_models/plotting_data/",
        fdr="{fdr}",
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
            "phenotype_prediction/models/linear_models/plotting_data/all_recomputed_metrics_test_{phenotype_suffix}_{fdr}.Rds",
            phenotype_suffix=phenotype_suffixes,
            fdr=fdrs,
        ),
        expand(
            "phenotype_prediction/models/logistic_models/plotting_data/combined_metrics_{phenotype_suffix}_{fdr}.Rds",
            phenotype_suffix=phenotype_suffixes,
            fdr=fdrs,
        ),


rule prep_metrics_linear_model:
    conda:
        "r-env"
    input:
        expand(
            "phenotype_prediction/models/linear_models/{phenotype}/{{phenotype_suffix}}_fdr-{{fdr}}.Rds",
            phenotype=phenotypes,
        ),
    output:
        "phenotype_prediction/models/linear_models/plotting_data/all_recomputed_metrics_test_{phenotype_suffix}_{fdr}.Rds",
    params:
        input_res_path="phenotype_prediction/models/linear_models/",
        out_dir="phenotype_prediction/models/linear_models/plotting_data/",
        fdr="{fdr}",
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
        expand(
            "phenotype_prediction/models/logistic_models/{phenotype}/{{phenotype_suffix}}_{top_bottom_q}-{top_q}_fdr-{{fdr}}.Rds",
            phenotype=phenotypes,
            top_q=top_quantiles,
            top_bottom_q=top_bottom_quantiles,
        ),
    output:
        metrics="phenotype_prediction/models/logistic_models/plotting_data/combined_metrics_{phenotype_suffix}_{fdr}.Rds",
    params:
        input_res_path="phenotype_prediction/models/logistic_models/",
        out_dir="phenotype_prediction/models/logistic_models/plotting_data/",
        fdr="{fdr}",
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
            "phenotype_prediction/models/linear_models/{phenotype}/{phenotype_suffix}_fdr-{fdr}.Rds",
            phenotype=phenotypes,
            phenotype_suffix=phenotype_suffixes,
            fdr=fdrs,
        ),
        expand(
            "phenotype_prediction/models/logistic_models/{phenotype}/{phenotype_suffix}_{top_bottom_q}-{top_q}_fdr-{fdr}.Rds",
            phenotype=phenotypes,
            phenotype_suffix=phenotype_suffixes,
            top_q=top_quantiles,
            fdr=fdrs,
            top_bottom_q=top_bottom_quantiles,
        ),


rule fit_logistic_regression_model_fdr:
    conda:
        "r-env"
    input:
        expand(
            "phenotype_prediction/r_data/{{phenotype}}/cv_split{cvsplit}/data.finished",
            cvsplit=range(cv_splits),
        ),
    output:
        "phenotype_prediction/models/logistic_models/{phenotype}/{phenotype_suffix}_{top_bottom_q}-{top_q}_fdr-{fdr}.Rds",
    resources:
        mem_mb=lambda wildcards, attempt: 20480 + (attempt - 1) * 4098,
        load=16000,
    params:
        model_out_file="phenotype_prediction/models/logistic_models/{phenotype}/{phenotype_suffix}_{top_bottom_q}-{top_q}_fdr-{fdr}.Rds",
        phenotype="{phenotype}",
        phenotype_suffix="{phenotype_suffix}",
        top_q="{top_q}",
        regression_data_input_path="phenotype_prediction/r_data/{phenotype}/",
        config=regression_config["r_config"],
        top_bottom_q="{top_bottom_q}",
    script:
        f"{code_dir}/run_logistic_regression.R"


rule fit_linear_regression_model_fdr:
    conda:
        "r-env"
    input:
        expand(
            "phenotype_prediction/r_data/{{phenotype}}/cv_split{cvsplit}/data.finished",
            cvsplit=range(cv_splits),
        ),
    output:
        "phenotype_prediction/models/linear_models/{phenotype}/{phenotype_suffix}_fdr-{fdr}.Rds",
    resources:
        mem_mb=lambda wildcards, attempt: 20480 + (attempt - 1) * 4098,
        load=16000,
    params:
        model_out_file="phenotype_prediction/models/linear_models/{phenotype}/{phenotype_suffix}_fdr-{fdr}.Rds",
        phenotype="{phenotype}",
        phenotype_suffix="{phenotype_suffix}",
        regression_data_input_path="phenotype_prediction/r_data/{phenotype}/",
        config=regression_config["r_config"],
    script:
        f"{code_dir}/run_linear_regression.R"

rule all_r_input_data:
    input:
        expand("phenotype_prediction/r_data/{phenotype}/cv_split{cvsplit}/data.finished",
               phenotype = phenotypes, cvsplit=range(cv_splits))
rule prep_data_for_r:
    input:
        burdens=expand(
            "phenotype_prediction/burdens/cv_split{cvsplit}/{split}_{btype}_burdens.pickle",
            split=["train", "test"],
            cvsplit=range(cv_splits),
            btype=burden_btypes,
        ),
        config="config_eval.yaml",
        dataset_train="cv_split{cvsplit}/deeprvat/{phenotype}/deeprvat/association_dataset.pkl",
        dataset_test="cv_split{cvsplit}/deeprvat/{phenotype}/deeprvat/association_dataset_test.pkl",
        discoveries="cv_split{cvsplit}/deeprvat/eval/all_significant.parquet",
    output:
        "phenotype_prediction/r_data/{phenotype}/cv_split{cvsplit}/data.finished",
    params:
        dataset_dir="cv_split{cvsplit}/deeprvat/{phenotype}/deeprvat/",
        burdens_dir="phenotype_prediction/burdens/cv_split{cvsplit}",
        out_dir="phenotype_prediction/r_data/{phenotype}/cv_split{cvsplit}",
        burden_types=[f"--burden-type {b} " for b in btypes],
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
                    " --dataset-dir {params.dataset_dir} "
                    " --burdens-dir {params.burdens_dir} "
                    " --out-dir {params.out_dir} "
                    " {input.config}"
                ),
                "touch {output}",
            ]
        )


rule all_burdens:
    input:
        expand(
            "phenotype_prediction/burdens/cv_split{cvsplit}/{split}_{btype}_burdens.pickle",
            split=["train", "test"],
            cvsplit=range(cv_splits),
            btype=burden_btypes,
        ),


burden_dir_mapper = {
    btype: "alternative_burdens" if btype != "deeprvat" else btype
    for btype in burden_btypes
}
split_mapper = {"train": "", "test": "_test"}
burden_phenotype = phenotypes[0]


# exctract burdens/ DeepRVAT gene impairment scores for significant genes
rule extract_burdens:
    input:
        discoveries="cv_split{cvsplit}/deeprvat/eval/all_significant.parquet",
        burdens=lambda wildcards: f"cv_split{{cvsplit}}/{burden_dir_mapper[wildcards.btype]}/{burden_phenotype}/{wildcards.btype}/burdens{split_mapper[wildcards.split]}/burdens.zarr",
    output:
        "phenotype_prediction/burdens/cv_split{cvsplit}/{split}_{btype}_burdens.pickle",
    resources:
        mem_mb=16000,
        load=16000,
    params:
        burden_input_dir=lambda wildcards: f"cv_split{wildcards.cvsplit}/{burden_dir_mapper[wildcards.btype]}/{burden_phenotype}/{wildcards.btype}/burdens{split_mapper[wildcards.split]}",
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


rule combine_discoveries:
    input:
        discoveries=expand(
            "cv_split{{cv_split}}/deeprvat/{phenotype}/deeprvat/eval/significant.parquet",
            phenotype=phenotypes,
        ),
    output:
        "cv_split{cv_split}/deeprvat/eval/all_significant.parquet",
    resources:
        mem_mb=8000,
        load=16000,
    params:
        sig_files=lambda wildcards, input: "".join(
            [f"--sig-file {d} " for d in input.discoveries]
        ),
        training_phenotype="".join(
            [f"--train-pheno {d} " for d in training_phenotypes]
        ),
    shell:
        " && ".join(
        [
            conda_check,
            (
        py + "extract_burdens_for_prs.py combine-significant "
                    " {params.sig_files} "
                    " {params.training_phenotype} "
                    " {output}"
                ),
            ]
        )


# include: '../cv_deeprvat_training/run_deeprvat_cv.snakefile'
