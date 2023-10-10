require(arrow)
require(plyr)
require(dplyr)
require(stringr)
require(tidyr)
require(data.table)

phenotype_suffix = snakemake@params[["phenotype_suffix"]] #TODO delete this later

code_dir = snakemake@params[["code_dir"]]
out_path = snakemake@output[["out_path"]]
query_phenotypes = snakemake@params[["phenotypes"]]
print(paste("Analyzing results for phenotypes:", query_phenotypes))


source(file.path(code_dir, "../../utils.R"))  #get phenotype renamer

phenotypes_map = gsub("_", " ", query_phenotypes)
names(phenotypes_map) = query_phenotypes
to_add <- query_phenotypes[!(query_phenotypes %in% names(phenotypes_deeprvat))]
phenotypes_deeprvat <- c(phenotypes_deeprvat, mapply(function(name, value) setNames(list(value), name),
                                       to_add, phenotypes_map[to_add]))

print(paste('phenotype map', phenotypes_deeprvat))


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
            print(paste(this_res_file, 'exists'))
            this_df = read_parquet(this_res_file)
            all_baseline_res = rbind(all_baseline_res, this_df)
        } else {
            # file_phenotype = paste0(str_remove(pheno, "_standardized"), "_standardized")
            this_res_file = file.path(base_exp_dir, file_phenotype, "eval", "postprocessed_associations_testing.parquet")
            print(this_res_file)
            if (file.exists(this_res_file)) {
                print(paste(this_res_file, 'exists'))
                this_df = read_parquet(this_res_file) %>% rename('EAC_filtered'='CAF_filtered')
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
print(test_results %>% distinct(Trait))
test_results = test_results %>%
    drop_na(pval)
test_results = test_results %>%
    separate(Method, into = c("vtype", "ttype"), sep = "-", remove = FALSE)
test_results %>%
    distinct(Method)

eac_threshold = 5


pheno_col = "Trait"
pval_thres = 0.05

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

print("Computing significant couunts")
all_sig = test_results %>%
    rename(gene = gene_id) %>%
    filter(EAC_filtered >= eac_threshold) %>%
    filter(pval_type == "score")
all_sig = all_sig %>%
      group_by(across(pheno_col), pval_type) %>%
      mutate(qval = p.adjust(pval, method = "BH")) %>%
      mutate(Significant = qval < pval_thres) %>%
      ungroup() 
all_sig_counts = all_sig %>%
  drop_na(Significant) %>% 
#   filter(Significant) %>%
  distinct(Trait, Significant, gene, .keep_all = TRUE) %>%
  group_by(Trait) %>% 
  mutate(Significant = as.logical(Significant)) %>%
  summarise(discoveries = sum(Significant))

all_sig_counts = rbind(all_sig_counts, tibble(Trait = 'All traits', discoveries = sum(all_sig_counts[['discoveries']])))

saveRDS(all_sig, file.path(dirname(out_path), paste0('monti_results', phenotype_suffix, '.Rds')))
write_parquet(all_sig, file.path(dirname(out_path), paste0('monti_results', phenotype_suffix, '.parquet')))
saveRDS(all_sig_counts, file.path(dirname(out_path), paste0('monti_counts', phenotype_suffix, '.Rds')))
write_parquet(all_sig_counts, file.path(dirname(out_path), paste0('monti_counts', phenotype_suffix, '.parquet')))
  
print("Checking replication")
missing_replication_traits = setdiff(c(unlist(unique(test_results[pheno_col]))), unlist(unique(comparison[pheno_col])))
print(paste('Traits not present in replication data', missing_replication_traits, 'excluding these traits'))
test_results = as.data.table(test_results)
test_results = test_results[!(Trait %in% missing_replication_traits)]
test_results = as_tibble(test_results)
print(nrow(test_results %>% distinct(Trait)))

replication_monti <- add_replication(test_results %>%
    rename(gene = gene_id) %>%
    filter(EAC_filtered >= eac_threshold) %>%
    filter(pval_type == "score"), comparison, thresh = pval_thres, pheno_col = pheno_col) %>%
    select(Rank, Replicated, Significant, Method, Trait, pval_type, gene, {{pheno_col}}) #%>%
    # mutate(Trait = replace(Trait, Trait == "Mean platelet thrombocyte volume", "MPTVS"))
replication_monti %>%
    distinct(Trait)


print("Saving output")

saveRDS(replication_monti, out_path)
write_parquet(replication_monti, paste(str_split(out_path, "\\.")[[1]][1], "parquet", sep = "."))

print("Finished")