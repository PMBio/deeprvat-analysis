import os

from pathlib import Path
import pandas as pd


configfile: "config.yaml"


debug_flag = config.get("debug", False)
debug = "--debug " if debug_flag else ""

BTYPES_DICT = {
    "is_plof": "plof",
    "alphamissense": "alphamissense",
    "CADD_raw": "cadd",
    "AbSplice_DNA": "absplice",
    "PrimateAI_score": "primateai",
    "sift_score": "sift",
    "polyphen_score": "polyphen",
    "SpliceAI_delta_score": "splicai",
    "Consequence_missense_variant": "missense",
}


PLOF_CONSEQUENCES = [
    f"Consequence_{c}"
    for c in (
        "splice_acceptor_variant",
        "splice_donor_variant",
        "frameshift_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
    )
]

conda_check = 'conda info | grep "active environment"'
cuda_visible_devices = "echo CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"


wildcard_constraints:
    repeat="\d+",


alt_burdens_chunks = 30


btypes = config["alt_burdens_data"]["dataset_config"]["rare_embedding"]["config"][
    "annotations"
]

btypes = list(set(btypes) - set(PLOF_CONSEQUENCES))

btypes = [BTYPES_DICT[btype] if btype in BTYPES_DICT else btype for btype in btypes]
btypes.append("plof")


rule all:
    input:
        burdens_test=expand(
            "{btype}/burdens/chunk{altchunk}.finished",
            btype=btypes,
            altchunk=range(alt_burdens_chunks),
        ),


rule compute_alt_burdens_burdens:
    priority: 10
    input:
        config="config.yaml",
        dataset="{btype}/association_dataset.pkl",
    output:
        "{btype}/burdens/chunk{altchunk}.finished",
    params:
        prefix=".",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: 64000 + (attempt - 1) * 2 * 8098,
        load=8000,
        # gpus = 1
    shell:
        " && ".join(
        [
            conda_check,
            cuda_visible_devices,
            (
        "python alternative_burdens.py compute-alternative-burdens "
                    + debug
                    + " --n-chunks "
                    + str(alt_burdens_chunks)
                    + " "
                    "--chunk {wildcards.altchunk} "
                    "--dataset-file {input.dataset} "
                    "--btype {wildcards.btype} "
                    "{input.config} "
                    "{wildcards.btype}/burdens"
                ),
                "touch {output}",
            ]
        )


rule alt_burdens_association_dataset:
    input:
        config="config.yaml",
    output:
        "{btype}/association_dataset.pkl",
    threads: 4
    priority: 12
    resources:
        mem_mb=lambda wildcards, attempt: 32000 * (attempt + 1),
        load=64000,
    shell:
        " && ".join(
        [
            conda_check,
            (
        "python alternative_burdens.py make-dataset "
                    "--data-key alt_burdens_data "
                    "--btype {wildcards.btype} " + debug + "{input.config} "
                    "{output}"
                ),
            ]
        )
