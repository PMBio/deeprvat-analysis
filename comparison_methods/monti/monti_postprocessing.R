require(arrow)
require(plyr)
require(dplyr)
require(stringr)
require(tidyr)

cct_vtypes = c("missense", "spliceai")
pval_cols = c(score = "pv_score")
tests_to_exclude = c("plof-skat", "deepripe-burden")  #plof_skat and deepripe_burden is not done by Monti



config = snakemake@params[["config"]]
phenotype = snakemake@params[["phenotype"]]
code_dir = snakemake@params[["code_dir"]]
gene_file = config[["data"]][["dataset_config"]][["gene_file"]]
input_files = snakemake@input[["testing_associations"]]
out_path = snakemake@output[["out_path"]]

source(file.path(code_dir, "monti_utils.R"))


print(gene_file)
genes <- read_parquet(gene_file) %>%
    select(!"__index_level_0__") %>%
    separate_wider_delim(gene, delim = ".", names = c("ensembl97_id", "id_suffix"))





ReadResultTable <- function(input_files, phenotype, cols_to_keep = c("gene", "EAC", "pval", "EAC_filtered")) {

    pheno_name = str_remove(phenotype, "_standardized")
    all_baseline_res = tibble()
    for (file in input_files) {
        vtype = str_split(file, "/", simplify = TRUE)[2]
        ttype = str_split(file, "/", simplify = TRUE)[3]
        this_res = read_parquet(file) %>%
            select(any_of(cols_to_keep)) %>%
            mutate(vtype = vtype, ttype = ttype, phenotype_std = phenotype, phenotype = pheno_name)
        all_baseline_res = rbind.fill(all_baseline_res, this_res)
    }

    all_baseline_res = all_baseline_res %>%
        unite(col = "Method", vtype, ttype, sep = "-") %>%
        rename(gene_id = "gene") %>%
        left_join(genes %>%
            select(id, ensembl97_id, gene_name), by = c(gene_id = "id"))
    return(all_baseline_res)
}


test_results = ReadResultTable(input_files = input_files, phenotype = phenotype) %>%
    rename(pv_score = "pval")

test_results_processed = tibble()
for (pval_name in names(pval_cols)) {
    print(paste("Using pvalues from", pval_name, "as pvals"))
    pval_col = pval_cols[[pval_name]]
    # set column that you want to use as pvals
    this_test_results = test_results %>%
        mutate(`:=`(pval, !!as.name(pval_col)))
    this_test_results[[pval_col]][this_test_results[[pval_col]] >= 1] = 0.99
    this_test_results[['pval']][this_test_results[['pval']] >= 1] = 0.99

    cct_combined_pvals = cct_combine_pvals(test_results = this_test_results, cct_vtypes = cct_vtypes)

    cct_combined_pvals = cct_combined_pvals %>%
        left_join(test_results %>%
            select(Method, gene_id, EAC_filtered))

    this_test_results_processed = this_test_results %>%
        select(all_of(colnames(cct_combined_pvals))) %>%
        filter(!grepl(paste(cct_vtypes, collapse = "|"), Method)) %>%
        filter(!Method %in% tests_to_exclude) %>%
        rbind(cct_combined_pvals) %>%
        mutate(pval_type = pval_name)

    test_results_processed = rbind(test_results_processed, this_test_results_processed)

}
print("Writing output")

print(out_path)
write_parquet(test_results_processed, out_path)
print("succesfully completed")

