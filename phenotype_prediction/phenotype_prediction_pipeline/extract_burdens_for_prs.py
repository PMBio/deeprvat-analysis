# from prs_comparison import get_burdens
import pickle
import pandas as pd
import sys
import os
import logging
from pathlib import Path
import click
import yaml
import zarr
import numpy as np
from prs_comparison import PHENOTYPE_MAPPING

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def get_burdens(burden_dir, query_genes, burden_type):
    all_burdens = zarr.open(Path(burden_dir) / "burdens.zarr")  # [:, :, repeat]
    genes = pd.Series(np.load(Path(burden_dir) / "genes.npy"))
    logger.info(f"All burdens shape: {all_burdens.shape}")

    n_missing_genes = len(query_genes) - len(set(genes).intersection(query_genes))
    query_genes = set(genes).intersection(query_genes)

    logger.info(
        f"Number of genes query genes (genebass hits) that "
        f"are not in DeepRVAT results and therefore removed: {n_missing_genes}"
    )

    gene_id_index_dict = {gene: genes[genes == gene].index[0] for gene in query_genes}
    if burden_type == "deeprvat":
        n_repeats = all_burdens.shape[2]
        out_burdens = {
            gene: {repeat: {} for repeat in range(n_repeats)}
            for gene in gene_id_index_dict.keys()
        }
        for gene, index in gene_id_index_dict.items():
            logger.info(f"Retreiving burdens for gene {gene}")
            for repeat in range(n_repeats):
                burdens_this_gene = all_burdens.oindex[:, index, repeat]

                out_burdens[gene][repeat] = burdens_this_gene
    else:
        assert len(all_burdens.shape) == 2
        out_burdens = {}
        for gene, index in gene_id_index_dict.items():
            logger.info(f"Retreiving burdens for gene {gene}")
            burdens_this_gene = all_burdens.oindex[:, index]
            out_burdens[gene] = burdens_this_gene

    return out_burdens


@click.group()
def cli():
    pass


@cli.command()
@click.option("--discoveries-file", type=click.Path(exists=True))
@click.option("--burdens-out", type=click.Path(), default=None)
@click.option("--burden-input-dir", type=click.Path(exists=True))
@click.option("--burden-type", "-s")
@click.option("--debug", is_flag=True)
def extract_burdens(
    burden_input_dir, burdens_out, burden_type, discoveries_file, debug
):
    logger.info(f"Reading discoveries file {discoveries_file}")
    discoveries = pd.read_parquet(discoveries_file)  # TODO change this to 0.4

    logger.info(f"Retrieving for all phenotypes {burden_input_dir}")
    logger.info(f'Considered phenotypes are {discoveries["Trait"].unique()}')

    logger.info("Extracting genes based on FDR threshold")

    query_genes = set(discoveries["gene"])
    logger.info(f"Number of query genes: {len(query_genes)}")

    logger.info(f"Total number of genes: {len(query_genes)}")

    if debug:
        query_genes = list(query_genes)[:2]

    logger.info("Retreiving burdens")
    burden_dict = get_burdens(burden_input_dir, query_genes, burden_type)

    with open(burdens_out, "wb") as f:
        pickle.dump(burden_dict, f)

    logger.info("succesfully completed")


@cli.command()
@click.option("--sig-file", type=click.Path(exists=True), multiple=True)
@click.option("--min-discoveries", type=int, default=1)
@click.argument("out-file", type=click.Path())
def combine_significant(sig_file, min_discoveries, out_file):
    methods_to_keep = ["DeepRVAT", "Burden/SKAT combined"]
    phenotypes = [file.split("/")[2].replace("_", " ") for file in sig_file]

    logger.info(f"Combining discoveries for phenotypes {phenotypes}")
    sig_file_dict = dict(zip(phenotypes, sig_file))
    all_discoveries = pd.concat([pd.read_parquet(file) for file in sig_file])

    existing_combinations = {(p, m) for p in phenotypes for m in methods_to_keep}
    disc_counts = (
        all_discoveries.query("Method in @methods_to_keep")
        .groupby(["phenotype", "Method"])
        .size()
        .reset_index()
        .rename(columns={0: "discoveries"})
    )
    logger.info(f"Significant discoveries {disc_counts}")
    logger.info(
        f"Requiring at least {min_discoveries} discoveries. Taking genes with lowest pval of not enough genes are significant"
    )

    suff_counts = disc_counts.query("discoveries >= @min_discoveries")

    existing_data_set = set(
        tuple(row) for row in suff_counts[["phenotype", "Method"]].to_numpy()
    )
    logger.info(f"Combinations with sufficient discoveries: {suff_counts}")
    missing_combinations = existing_combinations - existing_data_set
    missing_combinations = pd.DataFrame(
        missing_combinations, columns=["phenotype", "Method"]
    )
    assert len(missing_combinations) + len(suff_counts) == len(existing_combinations)

    if min_discoveries == 0:
        logger.info("Min discoveries is set to 0. Writing significant genes only.")
        logger.info(
            "Caution: The following phenotype/method combinations have zero discoveries"
        )
        logger.info(missing_combinations)
    else:
        all_discoveries = suff_counts[["phenotype", "Method"]].merge(
            all_discoveries, how="left"
        )
        assert len(all_discoveries[["phenotype", "Method"]].drop_duplicates()) == len(
            suff_counts
        )

        for pheno in missing_combinations["phenotype"].unique():
            logger.info(f"Adding missing discoveries for phenotype {pheno}")
            missing_methods = missing_combinations.query("phenotype == @pheno")[
                "Method"
            ].unique()
            all_pvals = pd.read_parquet(
                f"{os.path.dirname(sig_file_dict[pheno])}/all_results.parquet"
            )
            top_associations = (
                all_pvals.query("Method in @missing_methods")
                .sort_values(["Method", "pval"])
                .drop_duplicates(subset=["phenotype", "gene", "Method"])
                .groupby("Method")
                .apply(lambda group: group.nsmallest(min_discoveries, "pval"))
                .reset_index(drop=True)
            )
            all_discoveries = pd.concat([all_discoveries, top_associations])

        assert len(all_discoveries[["phenotype", "Method"]].drop_duplicates()) == len(
            existing_combinations
        )

    all_discoveries = all_discoveries.rename(columns={"phenotype": "Trait"})
    this_pheno_mapping = dict(
        zip(
            [i.replace("_", " ") for i in PHENOTYPE_MAPPING.keys()],
            PHENOTYPE_MAPPING.values(),
        )
    )
    all_discoveries["Trait"] = all_discoveries["Trait"].replace(this_pheno_mapping)
    logger.info(
        f"Number of discoveries:  {all_discoveries.groupby(['Trait', 'Method']).size()}"
    )

    all_discoveries.to_parquet(out_file)


if __name__ == "__main__":
    cli()
