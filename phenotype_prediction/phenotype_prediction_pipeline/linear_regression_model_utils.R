
library(Rcpp)

library(stringr)
library(arrow)
library(yardstick)
library(ggplot2)
library(gdata)
library(dplyr)
library(tidyr)
library(logger)



FitLinearModel <- function(data_list, all_x_cols, model_name = "") {
    log_info("Fitting model")
    model_formula = paste0("y ~ ", paste(all_x_cols, collapse = " + "))
    print(paste0("Model formula: ", model_formula))
    model_data = data_list[["train"]][["model_data"]]
    model_data = cbind(model_data, as_tibble_col(unlist(data_list[["train"]][["Y"]]), column_name = "y")) %>%
        select(-y_bin)
    
    model = lm(model_formula, data = model_data)
    summary(model)

    # combined_res_train = EvalModel(model$fitted.values, unlist(data_list[["train"]][["Y_bin"]]), unlist(data_list[["train"]][["Y"]]), model_name, logistic_model = FALSE, split = "train")
    test_split = "test"
    fitted_results = predict(model, newdata = data_list[[test_split]][["model_data"]], type = "response")
    print(summary(model))
    combined_res = EvalModel(fitted_results, as.integer(data_list[[test_split]][["Y_bin"]]), data_list[[test_split]][["Y"]] %>%
        pull(var = 1), model_name, logistic_model = FALSE)

    # combined_res = FlattenDfList(list(combined_res_train, combined_res_test))
    combined_res = mapply(cbind, combined_res, model_name = model_name, SIMPLIFY = F)

    return_data = list(res = combined_res, model = model$coefficients)
    # return(combined_res)
    return(return_data)

}

EvalModel <- function(fitted_results, truth, y_true_cont, plot_title = "", decision_threshold = 0.5, logistic_model = TRUE) {
    fitted_results_bin <- ifelse(fitted_results > decision_threshold, 1, 0)
    this_res = tibble(estimate = fitted_results, estimate_bin = as.factor(fitted_results_bin), truth = as.factor(truth), Y = y_true_cont)
    if (logistic_model) {
        log_info("Computing metrics for Linear model")

        class_and_probs_metrics <- metric_set(roc_auc, pr_auc, accuracy, f_meas)
        this_metrics = this_res %>%
            class_and_probs_metrics(truth, estimate, estimate = estimate_bin, event_level = "second")
    } else {
        print("Computing metrics for linear model")
        this_metrics = metrics(this_res, Y, estimate)
    }

    res_list = list(res = this_res, metrics = this_metrics)
    return(res_list)
}


PredictPhenoRareBurdenLinear <- function(y_data, covariate_cols, phenotype_col, this_out_dir, top_q, gene_lists, btypes = c("plof", "missense"), genes_to_keep = NULL) {
    log_info("Fitting Phenotype model for Rare burden")
    all_res_list = list()
    all_models_list = list()
    for (btype in btypes) {
        print(btype)
        all_models_list[[btype]] = list()
        inner_list = list()
        for (gene_list in gene_lists) {
            # this_genes_to_keep = ifelse(is.null(genes_to_keep), NULL, list(genes_to_keep[[gene_list]]))
            if (is.null(genes_to_keep)) {
                this_genes_to_keep = NULL
            } else {
                this_genes_to_keep = list(genes_to_keep[[gene_list]])
            }
            all_model_data = LoadDataForModel(y_data, covariate_cols = covariate_cols, phenotype_col = phenotype_col, top_q = top_q, this_out_dir = this_out_dir, btype = btype, gene_list = gene_list, genes_to_keep = this_genes_to_keep)
            data_list = all_model_data[["data_list"]]
            all_x_cols = all_model_data[["all_x_cols"]]
            model_output = FitLinearModel(data_list, all_x_cols, model_name = paste(btype, gene_list, sep = "-"))

            all_res_list = append(all_res_list, list(model_output[["res"]]))
            all_models_list[[btype]][[gene_list]] = model_output[["model"]]

        }
    }
    all_res_list = FlattenDfList(all_res_list)
    return(all_res_list)
}

PredictPhenoCovariatesLinear <- function(y_data, covariate_cols, phenotype_col, this_out_dir, top_q, use_top_q = TRUE) {
    log_info("Fitting Phenotype model for covariates only")
    # use any btype here
    all_model_data = LoadDataForModel(y_data, covariate_cols = covariate_cols, phenotype_col = phenotype_col, top_q = top_q, this_out_dir = this_out_dir, btype = "deeprvat", gene_list = "baseline_only", use_top_q = use_top_q)  #can basically use any combination of method, key_1, key_2 since the covariates are all the same

    data_list = all_model_data[["data_list"]]
    all_x_cols = all_model_data[["all_x_cols"]]
    covariates_to_keep = c(age, sex, genetic_pcs)
    cols_to_keep = unlist(to_list(for (string in all_x_cols) if (any(sapply(covariates_to_keep, grepl, string)))
        string))

    all_res_list = list()
    for (useprs in c(FALSE, TRUE)) {
        # print(paste0('use prs ', useprs))
        this_cols_to_keep = cols_to_keep
        if (useprs) {
            if (!("prs" %in% all_x_cols)) {
                break
            } else {
                model_name = "covariates-PRS"
                this_cols_to_keep = append(this_cols_to_keep, "prs")
            }

        } else {
            model_name = "covariates"
        }
        # print(this_cols_to_keep)
        model_output = FitLinearModel(data_list, this_cols_to_keep, model_name = model_name)

        all_res_list = append(all_res_list, list(model_output[["res"]]))

        # names(inner_list) = gene_lists all_res_list = append(all_res_list, list(inner_list))
    }
    all_res_list = FlattenDfList(all_res_list)
    return(all_res_list)
}

irnt <- function(x) {
    numPhenos <- sum(!is.na(x))
    quantilePheno <- (rank(x) - 0.5)/numPhenos
    phenoIRNT <- qnorm(quantilePheno)
    return(phenoIRNT)
}
