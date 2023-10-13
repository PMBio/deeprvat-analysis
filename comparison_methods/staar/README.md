# Run STAAR 

Here we run STAAR (Li et al., Nature Genetics, 2020, https://github.com/xihaoli/STAAR) with the same continuous variant annotations as used for DeepRVAT. 

As suggested [here](https://github.com/xihaoli/STAARpipeline-Tutorial), we use five genetic categories to aggregate coding rare variants of each protein-coding gene: (1) putative loss of function (stop gain, stop loss and splice) RVs, (2) missense RVs, (3) disruptive missense RVs, (4) putative loss of function and disruptive missense RVs, and (5) synonymous RV.
P-values for the different masks are combined afterwards and corrected for multiple testing. 

To run the tests, generate an experiment directory with the `config.yaml`. 

## Input data
The experiment directory in addition requires to have the same input data as specified for DeepRVAT (https://github.com/PMBio/deeprvat/), including
- `annotations.parquet`
- `genes.parquet`
- `genotypes.h5`
- `variants.parquet`
- `phenotypes.parquet`


### Varaint annotations
As for the DeepRVAT input, the  `annotations.parquet` should have the variant id (`id`) as the index, and the columns `gene_ids` and `combined_UKB_NFE_AF`

Also, the PHRED score for all continuous annotations used with DeepRVAT have to be provided (`−10 × log10(rank(−score_a)/M)`, where M is total number of variants sequenced across the whole genome and the ranking is done for each annotation `a` seperately)
- CADD_PHRED
- sift_score_PHRED
- polyphen_score_PHRED
- condel_score_PHRED
- DeepSEA_PC_1_PHRED
- DeepSEA_PC_2_PHRED
- DeepSEA_PC_3_PHRED
- DeepSEA_PC_4_PHRED
- DeepSEA_PC_5_PHRED
- DeepSEA_PC_6_PHRED
- PrimateAI_score_PHRED
- AbSplice_DNA_PHRED
- DeepRipe_plus_QKI_lip_hg2_PHRED
- DeepRipe_plus_QKI_clip_k5_PHRED
- DeepRipe_plus_KHDRBS1_clip_k5_PHRED
- DeepRipe_plus_ELAVL1_parclip_PHRED
- DeepRipe_plus_TARDBP_parclip_PHRED
- DeepRipe_plus_HNRNPD_parclip_PHRED
- DeepRipe_plus_MBNL1_parclip_PHRED
- DeepRipe_plus_QKI_parclip_PHRED
- SpliceAI_delta_score_PHRED   

Also, to build the variant filter masks, the following columns have to be provided. 

| Column name | Description  |
|:----------|:----------|
| is_plof    | Binary indicator. 1 if any VEP consquence out of  splice_acceptor_variant, splice_donor_variant, frameshift_variant, stop_gained, stop_lost, start_lost is True  |
| Consequence_missense_variant    | VEP annoation    |
| disruptive_missense    | variant is annotated with `sift_prediction_deleterious` and `polyphen_prediction_probably_damaging` by VEP   | 
| plof_or_disruptive_missense    | `disruptive_missense` or `is_plof`    |
| Consequence_synonymous_variant    | VEP annoation    |


## Running the pipeline
From the experiment directories (`experiments/staar`/`experiments/binary`), then run `snakemake {this_directoy_path}/staar.snakefile --use-conda' and any other flags you might need.
The pipeline requires an R environment (see `r-env.yaml` file). Also, STAAR has to be installed 
`library(devtools)`
`devtools::install_github("xihaoli/STAAR")`

If the R environment name is different from `r-env`, The environment name has to be changed accordingly in `staar.snakefile` (currently `r-env`)
