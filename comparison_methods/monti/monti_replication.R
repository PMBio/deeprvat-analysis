require(arrow)
require(plyr)
require(dplyr)
require(stringr)
require(tidyr)
require(data.table)


code_dir = snakemake@params[["code_dir"]]
out_path = snakemake@output[["out_path"]]
query_phenotypes = snakemake@params[["phenotypes"]]
print(paste("Analyzing results for phenotypes:", query_phenotypes))


source(file.path(code_dir, "../../utils.R"))  #get phenotype renamer
genes <- read_parquet("genes.parquet") %>%
    select(!"__index_level_0__") %>%
    separate_wider_delim(gene, delim = ".", names = c("ensembl97_id", "id_suffix"))


comparison <- read_parquet(file.path(code_dir, "../../data/comparison_results.parquet")) %>%
    mutate(in_comparison = TRUE) %>%
    mutate(Trait = str_remove(phenotype, "_standardized")) %>%
    mutate(Trait = recode(Trait, !!!phenotypes_deeprvat)) %>%
    rename(phenotype_standardized = phenotype) %>%
    mutate(phenotype = gsub("_standardized", "", phenotype_standardized))

ReadResultTables2 <- function(phenotypes) {
    all_baseline_res = tibble()
    base_exp_dir = "."
    for (pheno in phenotypes) {
        pheno_name = ifelse(pheno == "IGF_1", "IGF-1", pheno)
        file_phenotype = str_remove(pheno, "_standardized")
        print("Trying to read results for non-standardized phenotype")
        this_res_file = file.path(base_exp_dir, file_phenotype, "eval", "postprocessed_associations.parquet")
        print(this_res_file)
        if (file.exists(this_res_file)) {
            this_df = read_parquet(this_res_file)
            all_baseline_res = rbind(all_baseline_res, this_df)
        } else {
            file_phenotype = paste0(str_remove(pheno, "_standardized"), "_standardized")
            this_res_file = file.path(base_exp_dir, file_phenotype, "eval", "postprocessed_associations.parquet")
            print(this_res_file)
            if (file.exists(this_res_file)) {
                this_df = read_parquet(this_res_file)
                all_baseline_res = rbind(all_baseline_res, this_df)
            } else {
                message(paste("Failed to read results for phenotype", pheno, sep = ":"))
            }
        }
    }

    return(all_baseline_res)
}



print("reading results")
test_results = ReadResultTables2(phenotypes = query_phenotypes) %>%
    mutate(Trait = recode(phenotype, !!!phenotypes_deeprvat))

test_results = test_results %>%
    drop_na(pval)
test_results = test_results %>%
    separate(Method, into = c("vtype", "ttype"), sep = "-", remove = FALSE)
test_results %>%
    distinct(Method)

eac_threshold = 5
df = test_results %>%
    rename(gene = gene_id) %>%
    filter(EAC_filtered >= eac_threshold) %>%
    filter(pval_type == "score")

pheno_col = "Trait"

add_replication <- function(df, comparison, thresh = 0.05, pheno_col = "phenotype") {
    df <- df %>%
        group_by(across(pheno_col), pval_type) %>%
        mutate(qval = p.adjust(pval, method = "BH")) %>%
        mutate(Significant = qval < thresh) %>%
        ungroup()
    rep_list = list()

    for (this_pval_type in unique(df[["pval_type"]])) {
        this_df = df %>%
            filter(pval_type == this_pval_type)
        df_all <- this_df %>%
            arrange(qval) %>%
            distinct(across(pheno_col), gene, .keep_all = TRUE) %>%
            left_join(comparison, by = c(pheno_col, "gene")) %>%
            replace_na(list(in_comparison = FALSE)) %>%
            mutate(Rank = row_number(), Replicated = cumsum(in_comparison)) %>%
            mutate(`:=`(!!pheno_col, "all_phenotypes"), Trait = "All traits") %>%
            head(1000)

        rep_list = append(rep_list, list(df_all))
        for (this_pheno in unique(this_df[[pheno_col]])) {

            pheno_rep <- this_df %>%
                filter(across(all_of(pheno_col)) == this_pheno) %>%
                arrange(qval) %>%
                distinct(across(pheno_col), gene, .keep_all = TRUE) %>%
                left_join(comparison, by = c(pheno_col, "gene")) %>%
                replace_na(list(in_comparison = FALSE)) %>%
                mutate(Rank = row_number(), Replicated = cumsum(in_comparison)) %>%
                head(1000)
            rep_list = append(rep_list, list(pheno_rep))
        }
    }
    df_combined = rbindlist(rep_list, use.names = TRUE)
    return(df_combined)
}

print("Checking replication")

replication_monti <- add_replication(test_results %>%
    rename(gene = gene_id) %>%
    filter(EAC_filtered >= eac_threshold) %>%
    filter(pval_type == "score"), comparison, pheno_col = "Trait") %>%
    select(Rank, Replicated, Significant, Method, Trait, pval_type, gene) %>%
    mutate(Trait = replace(Trait, Trait == "Mean platelet thrombocyte volume", "MPTVS"))
replication_monti %>%
    distinct(Trait)


print("Saving output")

saveRDS(replication_monti, out_path)
write_parquet(replication_monti, paste(str_split(out_path, "\\.")[[1]][1], "parquet", sep = "."))

print("Finished")
