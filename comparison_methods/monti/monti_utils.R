cct_combine_pvals = function(test_results, cct_vtypes = c("missense", "spliceai")) {
    cct_combined_pvals = tibble()
    for (vtype in cct_vtypes) {
        this_cols = grep(paste0(vtype, "-"), names(test_results), value = TRUE)
        # use score pvalue to checkt if the p-value from the combined vartype + plof test should be used only if score-test pvalue for skat or burden <
        # 0.1
        do_cct_df = test_results %>%
            select(-pval, -EAC, -EAC_filtered) %>%
            pivot_wider(values_from = pv_score, names_from = Method)

        all_pval_cols = grep(vtype, colnames(do_cct_df), value = TRUE)
        cols_to_keep = grep("burden|skat", colnames(do_cct_df), value = TRUE, invert = TRUE)

        do_cct_df = do_cct_df %>%
            select(all_of(c(cols_to_keep, all_pval_cols))) %>%
            mutate(do_cct = ifelse(rowSums(select(., starts_with(paste0(vtype, "-"))) < 0.1) > 0, TRUE, FALSE)) %>%
            mutate(do_cct = replace_na(do_cct, FALSE))

        # debugging print('subsetting') test = test[1:100,]
        test = test_results %>%
            select(-starts_with("pv_"), -EAC, -EAC_filtered) %>%
            pivot_wider(values_from = pval, names_from = Method) %>%
            select(all_of(c(cols_to_keep, all_pval_cols))) %>%
            left_join(do_cct_df %>%
                select(-ends_with("burden"), -ends_with("skat")))
        ## cct combine for each test individually
        for (ttype in c("burden", "skat")) {
            # for (ttype in c('skat')){
            cols_for_cct = paste(c(vtype, paste(vtype, "plof", sep = "_")), ttype, sep = "-")
            print(cols_for_cct)

            cct_row <- function(row) {
                pvals <- row[1:(length(row) - 1)]
                do_cct = tail(row, 1)
                if (do_cct) {
                  if (all(complete.cases(pvals))) {
                    return(CCT(pvals))
                  } else {
                    return(NA)
                  }
                } else {
                  return(row[1])
                }
            }

            cct_value_name = paste("cct", vtype, ttype, sep = "-")
            test = test %>%
                rowwise() %>%
                mutate(`:=`(!!cct_value_name, cct_row(c_across(cols = all_of(c(cols_for_cct, "do_cct"))))))

        }
        print(paste0("Number of rows (should be = #genes):", nrow(test)))
        test = test %>%
            pivot_longer(cols = matches("burden|skat"), names_to = "Method", values_to = "pval") %>%
            filter(grepl("cct", Method))

        cct_combined_pvals = rbind(cct_combined_pvals, test)

    }
    cct_combined_pvals = cct_combined_pvals %>%
        select(-do_cct) %>%
        mutate(Method = str_remove(Method, "cct-"))
    return(cct_combined_pvals)
}

ReadResultTables <- function(exp_names, vtypes = c("missense", "plof"), ttypes = c("burden", "skat"), phenotype_dirs = phenotype_dirs_sub, cols_to_keep = c("gene",
    "EAC", "pval", "EAC_filtered")) {
    res_list = c()
    for (exp_name in names(exp_names)) {
        print(exp_name)
        base_exp_dir = exp_names[[exp_name]]
        for (pheno in phenotype_dirs) {
            pheno_name = ifelse(pheno == "IGF_1", "IGF-1", pheno)
            for (vtype in vtypes) {
                for (ttype in ttypes) {
                  tryCatch({
                    file_phenotype = str_remove(pheno, "_standardized")
                    print("Trying to read results for non-standardized phenotype")
                    this_res_dir = file.path(base_exp_dir, file_phenotype, vtype, ttype, "results", "burden_associations_testing.parquet")
                    this_df = read_parquet(this_res_dir) %>%
                      select(any_of(cols_to_keep)) %>%
                      mutate(vtype = vtype, ttype = ttype, phenotype = pheno_name, exp_name = exp_name)
                    res_list = append(res_list, list(this_df))

                  }, error = function(e1) {
                    tryCatch({
                      file_phenotype = paste0(str_remove(pheno, "_standardized"), "_standardized")
                      print("Trying to read results for standardized phenotype")
                      this_res_dir = file.path(base_exp_dir, file_phenotype, vtype, ttype, "results", "burden_associations_testing.parquet")
                      this_df = read_parquet(this_res_dir) %>%
                        select(any_of(cols_to_keep)) %>%
                        mutate(vtype = vtype, ttype = ttype, phenotype = pheno_name, exp_name = exp_name)
                      res_list = append(res_list, list(this_df))
                    }, error = function(e2) {
                      message(paste("Failed to read results for phenotype", pheno, sep = ":"))
                    })
                  })

                }
            }
        }
    }
    all_baseline_res = rbindlist(res_list, fill = TRUE)
    all_baseline_res = all_baseline_res %>%
        unite(col = "Method", vtype, ttype, sep = "-") %>%
        # mutate(method_old = Method) %>% mutate(Method = recode(Method, !!!method_recode)) %>%
    rename(gene_id = "gene") %>%
        left_join(genes %>%
            select(id, ensembl97_id, gene_name), by = c(gene_id = "id"))
    return(all_baseline_res)
}
