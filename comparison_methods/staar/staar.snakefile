from pathlib import Path
import os

configfile: "config.yaml"


debug_flag = config.get("debug", False)
phenotypes = config["phenotypes"]

n_chunks = config.get("n_chunks", 30) if not debug_flag else 2

debug = "--debug " if debug_flag else ""
persist_burdens = "--persist-burdens" if config.get("persist_burdens", False) else ""

conda_check = 'conda info | grep "active environment"'
# cuda_visible_devices = 'echo CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES'

DEEPRVAT_ANALYSIS_DIR = os.environ['DEEPRVAT_ANALYSIS_DIR']
pipeline_dir = f'{DEEPRVAT_ANALYSIS_DIR}/comparison_methods/staar/'
py_pipeline = f'python {pipeline_dir}'
r = f"Rscript {pipeline_dir}"


wildcard_constraints:
    repeat="\d+",
    chunk="\d+",


masks = [
    "plof",
    "missense",
    "disruptive_missense",
    "plof_disruptive_missense",
    "synonymous",
]

rule all:
    input:
        'replication_staar.Rds'

rule replication:
    conda:
        "r-env"
    input:
        expand(
            "{phenotype}/{mask}/results/burden_associations.parquet",
            phenotype=phenotypes,
            mask=masks,
        ),
    output:
        out_path="replication_staar.Rds",
    params:
        code_dir=pipeline_dir,
        phenotypes=phenotypes,
        masks=masks,
    threads: 1
    resources:
        mem_mb=16000,
        load=16000,
    script:
        f"{pipeline_dir}/staar_replication.R"



rule combine_regression_chunks:
    conda:
        "r-env"
    input:
        expand(
            "{{phenotype}}/{{mask}}/results/burden_associations_chunk{chunk}.Rds",
            chunk=range(n_chunks),
        ),
    output:
        "{phenotype}/{mask}/results/burden_associations.parquet",
    threads: 4
    resources:
        mem_mb=2048,
        load=2000,
    shell:
        " && ".join(
            [conda_check, r + "staar_combine_chunks.R  {input} {output}"]
        )


rule all_regression_results:
    priority: 100
    input:
        expand(
            "{phenotype}/{mask}/results/burden_associations_chunk{chunk}.Rds",
            phenotype=phenotypes,
            mask=masks,
            chunk=range(n_chunks),
        ),


rule regress:
    priority: 100
    conda:
        "r-env"
    input:
        genotype_file="{phenotype}/{mask}/data/genotypes_chunk{chunk}.h5",
        annotation_file="{phenotype}/{mask}/data/annotations_chunk{chunk}.h5",
        x_file="{phenotype}/{mask}/data/X_chunk{chunk}.h5",
        y_file="{phenotype}/{mask}/data/y_chunk{chunk}.h5",
    output:
        temp("{phenotype}/{mask}/results/burden_associations_chunk{chunk}.Rds"),
    threads: 4
    resources:
        # mem_mb = 16000,
        mem_mb=lambda wildcards, attempt: 16000 * 1 * (attempt + 1),
        load=32000,
    shell:
        " && ".join(
            [
                conda_check,
                r
                + "staar.R {input.genotype_file} {input.annotation_file} {input.x_file} {input.y_file} {output}",
            ]
        )


rule all_data:
    priority: 30
    input:
        # Each h5 file has one dataset per gene in the chunk
        genotype_file=expand(
            "{phenotype}/{mask}/data/genotypes_chunk{chunk}.h5",
            phenotype=phenotypes,
            mask=masks,
            chunk=range(n_chunks),
        ),
        annotation_file=expand(
            "{phenotype}/{mask}/data/annotations_chunk{chunk}.h5",
            phenotype=phenotypes,
            mask=masks,
            chunk=range(n_chunks),
        ),
        x_file=expand(
            "{phenotype}/{mask}/data/X_chunk{chunk}.h5",
            phenotype=phenotypes,
            mask=masks,
            chunk=range(n_chunks),
        ),
        y_file=expand(
            "{phenotype}/{mask}/data/y_chunk{chunk}.h5",
            phenotype=phenotypes,
            mask=masks,
            chunk=range(n_chunks),
        ),


rule build_data:
    input:
        data="{phenotype}/{mask}/association_dataset_full.pkl",
        dataset="{phenotype}/{mask}/association_dataset_pickled.pkl",
        config="{phenotype}/{mask}/config.yaml",
    output:
        genotype_file=temp("{phenotype}/{mask}/data/genotypes_chunk{chunk}.h5"),
        annotation_file=temp("{phenotype}/{mask}/data/annotations_chunk{chunk}.h5"),
        x_file=temp("{phenotype}/{mask}/data/X_chunk{chunk}.h5"),
        y_file=temp("{phenotype}/{mask}/data/y_chunk{chunk}.h5"),
    threads: 4
    priority: 30
    resources:
        mem_mb=32000,
        disk_mb=8000,
        load=64000,
    shell:
        " && ".join(
        [
            conda_check,
            (
                py
        + "staar_data.py build-data "
                    + debug
                    + " --n-chunks "
                    + str(n_chunks)
                    + " "
                    "--chunk {wildcards.chunk} "
                    "--dataset-file {input.dataset} "
                    "--data-file {input.data} "
                    " {input.config} "
                    "{output.genotype_file} "
                    "{output.annotation_file} "
                    "{output.x_file} "
                    "{output.y_file} "
                ),
            ]
        )


rule all_association_dataset:
    priority: 40
    input:
        expand(
            "{phenotype}/{mask}/association_dataset_full.pkl",
            phenotype=phenotypes,
            mask=masks,
        ),
        expand(
            "{phenotype}/{mask}/association_dataset_pickled.pkl",
            phenotype=phenotypes,
            mask=masks,
        ),


rule association_dataset:
    input:
        "{phenotype}/{mask}/config.yaml",
    output:
        full="{phenotype}/{mask}/association_dataset_full.pkl",
        pickled="{phenotype}/{mask}/association_dataset_pickled.pkl",
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
                py
        + "staar_data.py make-dataset "
                    + debug
                    + "--pickled-dataset-file {output.pickled} "  # outputs/dataset
                    "{input} "
                    "{output.full}"
                ),
            ]
        )


rule config:
    priority: 60
    input:
        config="config.yaml",
    output:
        "{phenotype}/{mask}/config.yaml",
    threads: 1
    resources:
        mem_mb=1024,
        load=1000,
    shell:
        " && ".join(
        [
            conda_check,
            (
                py
        + "staar_data.py update-config "
                    + "--phenotype {wildcards.phenotype} "
                    + "--variant-type {wildcards.mask} "
                    + " {input.config} "
                    + "{output}"
                ),
            ]
        )
