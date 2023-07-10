# Run kernel and collapsing tests as in Monti et al. 

Here, we re-implement the Kernel and Collapsing tests performed by Monti et al. (Nature Comm. 2023).

To run the tests, generate an experiment directory with the `config.yaml`. 

## Input data
The experiment directory in addition requires to have the same input data as specified for DeepRVAT (https://github.com/PMBio/deeprvat/), including
- `annotations.parquet`
- `genes.parquet`
- `genotypes.h5`
- `variants.parquet`
- `phenotypes.parquet`


The `annotations.parquet` should have the follwowing columns:

| Column name | Description  |
|:----------|:----------|
| id    | **Has to be the data frame index!** variant id as used  in `annotations.parquet` and `variants.parquet` |
| gene_ids    | gene id of the gene the variant is assigned to (as in `genes.parquet`|
| chrom    | chromosome   |
| pos   | variant position   |
| combined_UKB_NFE_AF   | variant MAF (Maximum of the MAF in the UK Biobank cohort and in gnomAD release 3.0 (non-Finnish European population)    |
| is_plof    | Binary indicator. 1 if any VEP consquence out of  splice_acceptor_variant, splice_donor_variant, frameshift_variant, stop_gained, stop_lost, start_lost is True  |
| Consequence_missense_variant    | VEP annoation    |
| missense_or_plof    | variant is annotated with `is_plof` or `Consequence_missense_variant` or both    | 
| polyphen_score    | Polyphen2 score    |
| sift_score    | SIFT score    |
| SpliceAI_delta_score    | Maximum of the four delta scores returned by SpliceAI   |
| keep_SpliceAI    | Binary indicator; 1 if   SpliceAI_delta_score > 0.1  |
| plof_or_keep_SpliceAI    | Binary indicator; 1 if   SpliceAI_delta_score > 0.1  or is_plof   |
| QKI_lip_hg2    | DeepRIPE prediction   |
| QKI_clip_k5    | DeepRIPE prediction     |
| TARDBP_parclip    |  DeepRIPE prediction   |
| HNRNPD_parclip    |  DeepRIPE prediction   |
| MBNL1_parclip    |  DeepRIPE prediction    |
| QKI_parclip    |  DeepRIPE prediction  |
| keep_DeepRipe    | Binary indicator; 1 if any DeepRIPE prediction has an abosolut value greater or equal to 0.25  |

### Variant to amino acid mapping
In addition, `amino_acid_variant_map.parquet` is required, which has the following columns
- `id` (should be the index column!)
- `gene_ids` single gene id
- `chrom`
- `pos`
- `Protein_position`: Specifies variant position in a protein (e.g., 3/100 denotes the variant at the 3rd amino acid out of 100). Any string/number can be used as long as variants in the same amino acid for a gene share a unique Protein_position. This is used to aggregate variants at the amino acid level

## Running the pipeline
From your experiemnt directory, then run `snakemake {this_directoy_path}/monti.snakefile --use-conda' and any other flags you might need.
The pipeline requires an r environment (see environment file). The environment name has to be changed accordingly in `monti.snakefile` (currently `r-env`)
