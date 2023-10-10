require(arrow)
require(plyr)
require(dplyr)
require(stringr)
require(tidyr)
require(ggplot2)
require(data.table)

phenotype_suffix = snakemake@params[["phenotype_suffix"]] #TODO delete this later

masks = snakemake@params[["masks"]]
code_dir = snakemake@params[["code_dir"]]
out_path = snakemake@output[["out_path"]]
query_phenotypes = snakemake@params[["phenotypes"]]
source(file.path(code_dir, "../../utils.R"))  #get phenotype renamer

print(paste("Analyzing results for phenotypes:", query_phenotypes))

phenotypes_map = gsub("_", " ", query_phenotypes)
names(phenotypes_map) = query_phenotypes
to_add <- query_phenotypes[!(query_phenotypes %in% names(phenotypes_deeprvat))]
phenotypes_deeprvat <- c(phenotypes_deeprvat, mapply(function(name, value) setNames(list(value), name),
                                       to_add, phenotypes_map[to_add]))
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
            print(paste('using this res file', this_res_file))
            this_df = read_parquet(this_res_file) %>%
                mutate(phenotype = file_phenotype, mask = mask)
            all_staar_res = rbind.fill(all_staar_res, this_df)
        } else {
            this_res_file = file.path(base_exp_dir, file_phenotype, mask, "results", "burden_associations_testing.parquet")
            print(paste('using this res file', this_res_file))
            this_df = read_parquet(this_res_file) %>%
            mutate(phenotype = file_phenotype, mask = mask) %>% rename('EAC'='CAF')
            all_staar_res = rbind.fill(all_staar_res, this_df)
        }
    }
}
print(all_staar_res %>% distinct(phenotype))
all_staar_res = all_staar_res %>%
    mutate(Trait = recode(phenotype, !!!phenotypes_deeprvat))

correction_method = "BH"
eac_threshold = 10  #same as used by Li et al. 2020 (STAAR paper)
staar_corrected = all_staar_res %>%
    filter(EAC >= eac_threshold) %>%
    group_by(Trait) %>%
    mutate(pval_adj = p.adjust(pval, method = correction_method)) %>%
    mutate(Significant = factor(pval_adj < alpha, levels = c("TRUE", "FALSE"))) %>%
    ungroup()

sig_counts_staar = staar_corrected %>%
    drop_na(Significant) %>% 
    distinct(Trait, Significant, gene, .keep_all = TRUE) %>%
    group_by(Trait) %>%
    mutate(Significant = as.logical(Significant)) %>%
    summarise(discoveries = sum(Significant))
    
sig_counts_staar = rbind(sig_counts_staar, 
    tibble(Trait = 'All traits',
    discoveries = sum(sig_counts_staar[['discoveries']])))

print(sig_counts_staar)
sum(sig_counts_staar %>%
    pull(discoveries))

saveRDS(staar_corrected, file.path(dirname(out_path), paste0('staar_corrected', phenotype_suffix, '.Rds')))
saveRDS(sig_counts_staar, file.path(dirname(out_path), paste0('sig_counts_staar', phenotype_suffix, '.Rds')))
write_parquet(sig_counts_staar, file.path(dirname(out_path), paste0('sig_counts_staar', phenotype_suffix, '.parquet')))

comparison <- read_parquet(file.path(code_dir, "../../data/comparison_results.parquet")) %>%
    mutate(in_comparison = TRUE) %>%
    mutate(Trait = str_remove(phenotype, "_standardized")) %>%
    mutate(Trait = recode(Trait, !!!phenotypes_deeprvat)) %>%
    rename(phenotype_standardized = phenotype) %>%
    mutate(phenotype = gsub("_standardized", "", phenotype_standardized))

add_replication <- function(df, comparison, thresh = 0.05, pheno_col = "phenotype") {
    df <- df %>%
        group_by(across(pheno_col)) %>%
        mutate(qval = p.adjust(pval, method = "BH")) %>%
        mutate(Significant = qval < thresh) %>%
        # mutate(Significant = factor(qval < thresh, levels = c('TRUE', 'FALSE'))) %>%
    ungroup()
    print("getting sig hits")
    rep_list = list()

    this_df = df
    print("all phenos")
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
    df_combined = rbindlist(rep_list, use.names = TRUE)
    return(df_combined)
}

# all_staar_res = all_staar_res %>%
#     left_join(comparison %>%
#         select(phenotype, Trait))
missing_replication_traits = setdiff(c(unlist(unique(all_staar_res['Trait']))), 
    unlist(unique(comparison['Trait'])))
print(paste('Traits not present in replication data', missing_replication_traits, 'excluding these traits'))
all_staar_res = as.data.table(all_staar_res)
all_staar_res = all_staar_res[!(Trait %in% missing_replication_traits)]
all_staar_res = as_tibble(all_staar_res)
print(all_staar_res %>% distinct(Trait))

# assertthat::assert_that(nrow(all_staar_res %>%
#     filter(is.na(Trait))) == 0)

if (nrow(all_staar_res %>%
    filter(is.na(Trait))) == 0){
        print('caution: not all Traits present in replicatoin data set!')
    }

replication_staar = add_replication(all_staar_res, comparison %>%
    select(-phenotype_standardized, phenotype), pheno_col = "Trait") %>%
    # select(-phenotype_standardized, phenotype, -`__index_level_0__`), pheno_col = "Trait") %>%
    rename(Method = mask) %>%
    select(Rank, Replicated, Significant, Method, Trait, gene) %>%
    mutate(exp_name = "STAAR")

print("Saving output")

saveRDS(replication_staar, out_path)
write_parquet(replication_staar, paste(str_split(out_path, "\\.")[[1]][1], "parquet", sep = "."))

print("Finished")
