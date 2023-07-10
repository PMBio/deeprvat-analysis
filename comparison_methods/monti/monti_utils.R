
#' An analytical p-value combination method using the Cauchy distribution
#'
#' The \code{CCT} function takes in a numeric vector of p-values, a numeric
#' vector of non-negative weights, and return the aggregated p-value using Cauchy method.
#' @param pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @param weights a numeric vector of non-negative weights. If \code{NULL}, the
#' equal weights are assumed (default = NULL).
#' @return the aggregated p-value combining p-values from the vector \code{pvals}.
#' @examples pvalues <- c(2e-02, 4e-04, 0.2, 0.1, 0.8)
#' @examples CCT(pvals = pvalues)
#' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
#' with analytic p-value calculation under arbitrary dependency structures.
#' \emph{Journal of the American Statistical Association}, \emph{115}(529), 393-402.
#' (\href{https://doi.org/10.1080/01621459.2018.1554485}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(3), 410-421.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.01.002}{pub})
#' @export
CCT <- function(pvals, weights = NULL) {
    ## copied from STAAR package check if there is NA
    if (sum(is.na(pvals)) > 0) {
        stop("Cannot have NAs in the p-values!")
    }

    #### check if all p-values are between 0 and 1
    if ((sum(pvals < 0) + sum(pvals > 1)) > 0) {
        stop("All p-values must be between 0 and 1!")
    }

    #### check if there are p-values that are either exactly 0 or 1.
    is.zero <- (sum(pvals == 0) >= 1)
    is.one <- (sum(pvals == 1) >= 1)
    if (is.zero && is.one) {
        stop("Cannot have both 0 and 1 p-values!")
    }
    if (is.zero) {
        return(0)
    }
    if (is.one) {
        warning("There are p-values that are exactly 1!")
        return(1)
    }

    #### check the validity of weights (default: equal weights) and standardize them.
    if (is.null(weights)) {
        weights <- rep(1/length(pvals), length(pvals))
    } else if (length(weights) != length(pvals)) {
        stop("The length of weights should be the same as that of the p-values!")
    } else if (sum(weights < 0) > 0) {
        stop("All the weights must be positive!")
    } else {
        weights <- weights/sum(weights)
    }

    #### check if there are very small non-zero p-values
    is.small <- (pvals < 1e-16)
    if (sum(is.small) == 0) {
        cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
    } else {
        cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
        cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
    }

    #### check if the test statistic is very large.
    if (cct.stat > 1e+15) {
        pval <- (1/cct.stat)/pi
    } else {
        pval <- 1 - pcauchy(cct.stat)
    }
    return(pval)
}


cct_combine_pvals = function(test_results, cct_vtypes = c("missense", "spliceai")) {
    cct_combined_pvals = tibble()
    for (vtype in cct_vtypes) {
        this_cols = grep(paste0(vtype, "-"), names(test_results), value = TRUE)
        # use score pvalue to checkt if the p-value from the combined vartype + plof test should be used only if score-test pvalue for skat or burden < 0.1
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

ReadResultTables <- function(exp_names, vtypes = c("missense", "plof"), ttypes = c("burden", "skat"), phenotype_dirs = phenotype_dirs_sub, cols_to_keep = c("gene", "EAC", "pval", "EAC_filtered")) {
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
