from pathlib import Path
import os

configfile: "config.yaml"


debug_flag = config.get("debug", False)

phenotypes = config["phenotypes"]

vtypes = config.get("variant_types", ["plof"])


ttypes = config.get("test_types", ["burden"])
rare_maf = config.get("rare_maf", 0.001)
n_chunks = config.get("n_chunks", 30) if not debug_flag else 2

debug = "--debug " if debug_flag else ""
persist_burdens = "--persist-burdens" if config.get("persist_burdens", False) else ""

DEEPRVAT_ANALYSIS_DIR = os.environ['DEEPRVAT_ANALYSIS_DIR']
pipeline_dir = f'{DEEPRVAT_ANALYSIS_DIR}/comparison_methods/monti/'
monti_pipeline = f"python {pipeline_dir}/monti_association.py "

conda_check = 'conda info | grep "active environment"'


wildcard_constraints:
    repeat="\d+",


rule all:
    input:
        'replication_monti.Rds'

rule replication:
    conda:
        "r-env"
    input:
        expand(
            "{phenotype}/eval/postprocessed_associations.parquet",
            phenotype=phenotypes,
        ),
    output:
        out_path="replication_monti.Rds",
    params:
        code_dir=pipeline_dir,
        phenotypes=phenotypes
    threads: 1
    resources:
        mem_mb=16000,
        load=16000,
    script:
        f"{pipeline_dir}/monti_replication.R"

rule postprocess_results:
    conda:
        "r-env"
    input:
        testing_associations=expand(
            "{{phenotype}}/{vtype}/{ttype}/results/burden_associations.parquet",
            vtype=vtypes,
            ttype=ttypes,
        ),
        config="{phenotype}/plof/config.yaml",  #simply any config
    output:
        out_path="{phenotype}/eval/postprocessed_associations.parquet",
    params:
        phenotype="{phenotype}",
        config=config,
        code_dir=pipeline_dir
    threads: 1
    resources:
        mem_mb=16000,
        load=16000,
    script:
        f"{pipeline_dir}/monti_postprocessing.R"


rule all_regression:
    input:
        expand(
            "{phenotype}/{vtype}/{ttype}/results/burden_associations.parquet",
            phenotype=phenotypes,
            vtype=vtypes,
            ttype=ttypes,
        ),


rule combine_regression_chunks:
    input:
        train=expand(
            "{{phenotype}}/{{vtype}}/{{ttype}}/results/burden_associations_chunk{chunk}.parquet",
            chunk=range(n_chunks),
        ),
    output:
        train="{phenotype}/{vtype}/{ttype}/results/burden_associations.parquet",
    threads: 1
    resources:
        mem_mb=2048,
        load=2000,
    shell:
        " && ".join(
            [
                conda_check,
                monti_pipeline + "combine-results " "{input.train} " "{output.train}",
            ]
        )


rule all_regression_results:
    input:
        expand(
            "{phenotype}/{vtype}/{ttype}/results/burden_associations_chunk{chunk}.parquet",
            phenotype=phenotypes,
            vtype=vtypes,
            ttype=ttypes,
            chunk=range(n_chunks),
        ),


rule regress:
    input:
        data="{phenotype}/{vtype}/association_dataset_full.pkl",
        dataset="{phenotype}/{vtype}/association_dataset_pickled.pkl",
        config="{phenotype}/{vtype}/config.yaml",
    output:
        out_path=temp(
            "{phenotype}/{vtype}/{ttype}/results/burden_associations_chunk{chunk}.parquet"
        ),
    threads: 10
    priority: 30
    resources:
        # mem_mb=24000,
        mem_mb=lambda wildcards, attempt: 32000 * (attempt + 1),
        load=8000,
    shell:
        " && ".join(
        [
            conda_check,
            (
                monti_pipeline
        + "run-association "
                    + debug
                    + " --n-chunks "
                    + str(n_chunks)
                    + " "
                    "--chunk {wildcards.chunk} "
                    "--dataset-file {input.dataset} "
                    "--data-file {input.data} " + persist_burdens + " "
                    " {input.config} "
                    "{wildcards.vtype} "
                    "{wildcards.ttype} "
                    "{output.out_path}"
                ),
            ]
        )


rule all_association_dataset:
    priority: 40
    input:
        expand(
            "{phenotype}/{vtype}/association_dataset_full.pkl",
            phenotype=phenotypes,
            vtype=vtypes,
        ),
        expand(
            "{phenotype}/{vtype}/association_dataset_pickled.pkl",
            phenotype=phenotypes,
            vtype=vtypes,
        ),


rule association_dataset:
    input:
        "{phenotype}/{vtype}/config.yaml",
    output:
        full="{phenotype}/{vtype}/association_dataset_full.pkl",
        pickled="{phenotype}/{vtype}/association_dataset_pickled.pkl",
    threads: 1
    priority: 40
    resources:
        mem_mb=40000,
        load=16000,
    shell:
        " && ".join(
        [
            conda_check,
            (
                monti_pipeline
        + "make-dataset "
                    + debug
                    + "--pickled-dataset-file {output.pickled} "
                    "{input} "
                    "{output.full}"
                ),
            ]
        )


rule config:
    input:
        config="config.yaml",
    output:
        "{phenotype}/{vtype}/config.yaml",
    params:
        rare_maf=str(rare_maf),
    threads: 1
    resources:
        mem_mb=1024,
        load=1000,
    shell:
        " && ".join(
        [
            conda_check,
            (
                monti_pipeline
        + "update-config "
                    + "--phenotype {wildcards.phenotype} "
                    + "--variant-type {wildcards.vtype} "
                    + "--rare-maf "
                    + "{params.rare_maf}"
                    + " {input.config} "
                    + "{output}"
                ),
            ]
        )
