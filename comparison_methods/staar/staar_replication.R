require(arrow)
require(plyr)
require(dplyr)
require(stringr)
require(tidyr)
require(ggplot2)
require(data.table)

phenotype_suffix = snakemake@params[["phenotype_suffix"]]  #TODO delete this later

masks = snakemake@params[["masks"]]
code_dir = snakemake@params[["code_dir"]]
out_path = snakemake@output[["out_path"]]
query_phenotypes = snakemake@params[["phenotypes"]]
correction_method = snakemake@params[["correction_method"]]
correction_method = ifelse(is.null(correction_method) | length(correction_method) == 0, "bonferroni", correction_method)

combine_pvals = snakemake@params[["combine_pvals"]]
combine_pvals = ifelse(length(combine_pvals) == 0 | is.null(combine_pvals), "bonferroni", combine_pvals)
print(sprintf("using p-value correction method %s", correction_method))

source(file.path(code_dir, "../../utils.R"))  #get phenotype renamer
source(file.path(code_dir, "../utils.R"))  #p-value combination functions 

# query_phenotypes = OLD_PHENOTYPES #TODO delete this, only for testing

print(paste("Analyzing results for phenotypes:", query_phenotypes))

phenotypes_map = gsub("_", " ", query_phenotypes)
names(phenotypes_map) = query_phenotypes
to_add <- query_phenotypes[!(query_phenotypes %in% names(phenotypes_deeprvat))]
phenotypes_deeprvat <- c(phenotypes_deeprvat, mapply(function(name, value) setNames(list(value), name), to_add, phenotypes_map[to_add]))
alpha = 0.05
all_staar_res = tibble()


for (pheno in query_phenotypes) {
    for (mask in masks) {
        pheno_name = ifelse(pheno == "IGF_1", "IGF-1", pheno)
        file_phenotype = str_remove(pheno, "_standardized")
        base_exp_dir = "."
        this_res_file = file.path(base_exp_dir, file_phenotype, mask, "results", "burden_associations.parquet")
        print(this_res_file)
        if (file.exists(this_res_file)) {
            print(paste("using this res file", this_res_file))
            this_df = read_parquet(this_res_file) %>%
                mutate(phenotype = file_phenotype, mask = mask)  # Added newly! seems like staar discoveries were in there twice! THIS MOST BE KEPT!!!!!!!
        } else {
            this_res_file = file.path(base_exp_dir, file_phenotype, mask, "results", "burden_associations_testing.parquet")
            print(paste("using this res file", this_res_file))
            this_df = read_parquet(this_res_file) %>%
                mutate(phenotype = file_phenotype, mask = mask) %>%
                rename(EAC = "CAF")
        }
        this_df = this_df %>%
            distinct() %>%
            mutate(pval = as.numeric(pval))
        all_staar_res = rbind.fill(all_staar_res, this_df)
    }
}
print(all_staar_res %>%
    distinct(phenotype))
eac_threshold = 10  #same as used by Li et al. 2020 (STAAR paper)
pheno_col = "Trait"

print(colnames(all_staar_res))
all_staar_res = all_staar_res %>%
    mutate(Trait = recode(phenotype, !!!phenotypes_deeprvat)) %>%
    filter(EAC >= eac_threshold)

print("aggregating pvals to gene level")
print(all_staar_res %>%
    group_by(gene, across(all_of(pheno_col))) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    arrange(n) %>%
    distinct(n))
if (!is.na(combine_pvals)) {
    all_staar_res = aggPvalsToGene(all_staar_res, combine_pvals, grouping_cols = c("gene", pheno_col)) %>%
        select(-pval_agg_method)
}
stopifnot(as.numeric(all_staar_res %>%
    group_by(gene, across(all_of(pheno_col))) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    arrange(n) %>%
    distinct(n)) == 1)

staar_corrected = all_staar_res %>%
    group_by(Trait) %>%
    mutate(p_adjust = p.adjust(pval, method = correction_method)) %>%
    mutate(Significant = factor(p_adjust < alpha, levels = c("TRUE", "FALSE"))) %>%
    ungroup()

sig_counts_staar = staar_corrected %>%
    drop_na(Significant) %>%
    distinct(Trait, Significant, gene, .keep_all = TRUE) %>%
    group_by(Trait) %>%
    mutate(Significant = as.logical(Significant)) %>%
    summarise(discoveries = sum(Significant))

sig_counts_staar = rbind(sig_counts_staar, tibble(Trait = "All traits", discoveries = sum(sig_counts_staar[["discoveries"]])))

print(sig_counts_staar)
sum(sig_counts_staar %>%
    filter(Trait == "All traits") %>%
    pull(discoveries))

saveRDS(staar_corrected, file.path(dirname(out_path), paste0("staar_corrected", phenotype_suffix, ".Rds")))
saveRDS(sig_counts_staar, file.path(dirname(out_path), paste0("sig_counts_staar", phenotype_suffix, ".Rds")))
write_parquet(sig_counts_staar, file.path(dirname(out_path), paste0("sig_counts_staar", phenotype_suffix, ".parquet")))

comparison <- read_parquet(file.path(code_dir, "../../data/comparison_results.parquet")) %>%
    mutate(in_comparison = TRUE) %>%
    mutate(Trait = str_remove(phenotype, "_standardized")) %>%
    mutate(Trait = recode(Trait, !!!phenotypes_deeprvat)) %>%
    rename(phenotype_standardized = phenotype) %>%
    mutate(phenotype = gsub("_standardized", "", phenotype_standardized))




add_replication <- function(df, comparison, thresh = 0.05, pheno_col = "phenotype") {
    df <- df %>%
        group_by(across(all_of(pheno_col))) %>%
        mutate(p_adjust = p.adjust(pval, method = correction_method)) %>%
        mutate(Significant = p_adjust < thresh) %>%
        # mutate(Significant = factor(p_adjust < thresh, levels = c('TRUE', 'FALSE'))) %>%
    ungroup()
    print("getting sig hits")
    rep_list = list()

    this_df = df
    print("all phenos")
    df_all <- this_df %>%
        arrange(pval) %>%
        distinct(across(all_of(pheno_col)), gene, .keep_all = TRUE) %>%
        left_join(comparison, by = c(pheno_col, "gene")) %>%
        replace_na(list(in_comparison = FALSE)) %>%
        mutate(Rank = row_number(), Replicated = cumsum(in_comparison)) %>%
        mutate(`:=`(!!pheno_col, "all_phenotypes"), Trait = "All traits") %>%
        head(1000)

    rep_list = append(rep_list, list(df_all))
    for (this_pheno in unique(this_df[[pheno_col]])) {
        pheno_rep <- this_df %>%
            filter(across(all_of(pheno_col)) == this_pheno) %>%
            arrange(pval) %>%
            distinct(across(pheno_col), gene, .keep_all = TRUE) %>%
            left_join(comparison, by = c(pheno_col, "gene")) %>%
            replace_na(list(in_comparison = FALSE)) %>%
            mutate(Rank = row_number(), Replicated = cumsum(in_comparison)) %>%
            head(1000)
        rep_list = append(rep_list, list(pheno_rep))

    }
    df_combined = rbindlist(rep_list, use.names = TRUE)
    return(df_combined)
}

# all_staar_res = all_staar_res %>% left_join(comparison %>% select(phenotype, Trait))
missing_replication_traits = setdiff(c(unlist(unique(all_staar_res["Trait"]))), unlist(unique(comparison["Trait"])))
warning(paste("Traits not present in replication data", missing_replication_traits, "excluding these traits"))


if (phenotype_suffix == "_training_phenotypes") {
    stopifnot(length(missing_replication_traits) == 0)
}

all_staar_res = as.data.table(all_staar_res)
all_staar_res = all_staar_res[!(Trait %in% missing_replication_traits)]
all_staar_res = as_tibble(all_staar_res)
print(all_staar_res %>%
    distinct(Trait))

# assertthat::assert_that(nrow(all_staar_res %>% filter(is.na(Trait))) == 0)


replication_staar = add_replication(all_staar_res, comparison %>%
    select(-phenotype_standardized, phenotype), pheno_col = "Trait") %>%
    mutate(Method = "STAAR") %>%
    select(Rank, Replicated, Significant, Method, Trait, gene) %>%
    # mutate(exp_name = 'STAAR') %>%
as_tibble()

print("Saving output")
out_file = str_split(out_path, "\\.")[[1]][1]

saveRDS(replication_staar, paste(out_file, "Rds", sep = "."))
write_parquet(replication_staar, paste(out_file, "parquet", sep = "."))


print("Finished")
