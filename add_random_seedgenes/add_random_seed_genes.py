import sys
import logging
from typing import Optional
import pandas as pd
import pickle
import yaml
import click

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level="INFO",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

@click.command()
@click.option("--exp_dir", type=click.Path(exists=True), default=".")
@click.option("--folds", type=int, default= 0, help="The number of different seed-gene setups, i.e. number of experiments each with a different set of added seed genes.")
@click.option("--p_genes_random", type=float, default= 0.2, help="Percentage of total seed genes to randomly add to each Fold.")
@click.option("--max_pval", type=float, default= 0.5, help="Select only seed genes to add with pvals <= max_pval.")
@click.option("--cv_splits", type=int, default= 5)
def seed_gene_add_selection(
    exp_dir: str,
    folds: int,
    p_genes_random: float, 
    max_pval: float,
    cv_splits: 5,
):
    '''
    This script is used for running the add-random-seed-genes analysis. 
    Works with the cv_split model setup.
    '''
    all_seed = {i : [] for i in range(folds)}
    added_gene_pval = {i : [] for i in range(folds)}

    for idx in range(folds):
        fold_dir = f'{exp_dir}/sg_set_{idx}'
        print(fold_dir)
        with open(f'{fold_dir}/deeprvat_config.yaml') as f:
            config = yaml.safe_load(f)
        phenotypes = config["training"]["phenotypes"]
        print(f"Adding Seed Genes to the following phenotypes:{phenotypes}")
        protein_coding_genes = pd.read_parquet(f'{fold_dir}/protein_coding_genes.parquet')

        for p in phenotypes: 
            print(p)
            seed_gene_file = f'{fold_dir}/base/cv_split0/deeprvat/{p}/deeprvat/seed_genes.parquet'
            this_seed = pd.read_parquet(seed_gene_file)
            if 'random_gene' in this_seed.columns:
                this_seed = this_seed.query('random_gene == False')
            else:
                this_seed = this_seed.assign(random_gene = False)
            n_genes_random = max(round(len(this_seed) * p_genes_random), 1)
            print(f'number of random genes that will be sampled: {n_genes_random}')

            # get baseline p-values of all genes
            this_base_file = f'{fold_dir}/base/baseline_results/{p}/eval/all_associations.parquet'
            this_base = pd.read_parquet(this_base_file)

            #sample new genes
            min_pval = this_base.query('EAC >= 50')[['gene', 'method', 'pval']].groupby('gene').min('pval').sort_values('pval').query(f'pval < {max_pval}')                                                                                                        
            new_genes = min_pval.reset_index()\
                .sample(n_genes_random)

            print(new_genes)

            gene_pval = new_genes.assign(phenotype=p)
            new_genes = new_genes['gene']
            seed_random = protein_coding_genes.query('id in @new_genes')\
                .assign(significant = True, correction_method = 'Bonferroni', random_gene = True)
                # addd columns such that they are the same as in this_seed
            seed_combined = pd.concat([this_seed, seed_random]).assign(phenotype = p)

            for split in range(cv_splits):
                seed_combined.to_parquet(f'{fold_dir}/cv_split{split}/deeprvat/{p}/deeprvat/seed_genes.parquet')

            print(len(seed_combined))
            all_seed[idx].append(seed_combined)
            added_gene_pval[idx].append(gene_pval)
        
        all_seed[idx] = pd.concat(all_seed[idx])
        added_gene_pval[idx] = pd.concat(added_gene_pval[idx])


    with open(f'{exp_dir}/added_gene_pval_dict.pkl','wb') as f:
        pickle.dump(added_gene_pval, f)

    with open(f'{exp_dir}/all_seed_dict.pkl','wb') as f:
        pickle.dump(all_seed, f)

if __name__ == "__main__":
    seed_gene_add_selection()