
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

    combined_res_train = EvalModel(model$fitted.values, unlist(data_list[["train"]][["Y_bin"]]), unlist(data_list[["train"]][["Y"]]), model_name, logistic_model = FALSE, split = "train")
    test_split = "test"
    fitted_results = predict(model, newdata = data_list[[test_split]][["model_data"]], type = "response")

    print(summary(model))
    combined_res_test = EvalModel(fitted_results, as.integer(data_list[[test_split]][["Y_bin"]]), data_list[[test_split]][["Y"]] %>%
        pull(var = 1), model_name, logistic_model = FALSE, split = "test")

    combined_res = FlattenDfList(list(combined_res_train, combined_res_test))
    combined_res = mapply(cbind, combined_res, model_name = model_name, SIMPLIFY = F)

    return_data = list(res = combined_res, model = model$coefficients)
    # return(combined_res)
    return(return_data)

}

EvalModel <- function(fitted_results, truth, y_true_cont, plot_title = "", decision_threshold = 0.5, logistic_model = TRUE, split = "test") {
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

    res_list = list(res = this_res %>%
        mutate(split = split), metrics = this_metrics %>%
        mutate(split = split))
    return(res_list)
}


PredictPhenoDeepRVATLinear <- function(y_data, covariate_cols, phenotype_col, this_out_dir, top_q, gene_lists, n_deeprvat_repeats = 6, genes_to_keep = NULL) {
    log_info("Fitting Phenotype model for DeepRVAT")
    all_res_list = list()
    all_models_list = list()
    for (gene_list in gene_lists) {
        log_info("Fitting model for gene_list: {gene_list}")
        all_models_list[[gene_list]] = list()
        estimate_list = list()

        if (is.null(genes_to_keep)) {
            this_genes_to_keep = NULL
        } else {
            this_genes_to_keep = list(genes_to_keep[[gene_list]])
        }

        for (deeprvat_repeat in seq(0, n_deeprvat_repeats - 1)) {
            print(paste("fitting model for deeprvat repeat", deeprvat_repeat))
            all_model_data = LoadDataForModel(y_data, covariate_cols = covariate_cols, phenotype_col = phenotype_col, top_q = top_q, this_out_dir = this_out_dir, method = "deeprvat", key_1 = gene_list, key_2 = deeprvat_repeat, genes_to_keep = this_genes_to_keep)
            data_list = all_model_data[["data_list"]]
            all_x_cols = all_model_data[["all_x_cols"]]
            model_output = FitLinearModel(data_list, all_x_cols, model_name = paste("deeprvat", deeprvat_repeat, gene_list, sep = "-"))

            all_res_list = append(all_res_list, list(model_output[["res"]]))

        }
    }

    all_res_list = FlattenDfList(all_res_list)
    return(all_res_list)
}

AverageDeepRVATBurdens <- function(all_data_list, key) {
    print(key)
    colnames = names(all_data_list[["repeat_0"]][["data_list"]][[key]][["model_data"]])
    stat_cols = grep("^gene_", colnames, value = TRUE, invert = TRUE)
    gene_cols = grep("^gene_", colnames, value = TRUE)

    first_df <- all_data_list[["repeat_0"]][["data_list"]][[key]][["model_data"]][gene_cols]
    # Initialize the average data frame with the same dimensions as the first data frame
    averages <- matrix(0, nrow = nrow(first_df), ncol = ncol(first_df))
    # Iterate over each 'repeat_x' value using a loop
    for (i in 0:5) {
        # Extract the data frame for the current 'repeat_x'
        current_df <- all_data_list[[paste0("repeat_", i)]][["data_list"]][[key]][["model_data"]][gene_cols]
        # Add the current data frame to the average data frame
        averages <- averages + current_df
    }
    # Calculate the average by dividing each cell by the number of data frames
    averages <- averages/6
    stat_cols = all_data_list[[paste0("repeat_", i)]][["data_list"]][[key]][["model_data"]][stat_cols]
    avg_df = cbind(stat_cols, averages)[colnames]

    return(avg_df)
}

PredictPhenoDeepRVATLinearAvgBurden <- function(y_data, covariate_cols, phenotype_col, this_out_dir, top_q, gene_lists, n_deeprvat_repeats = 6, genes_to_keep = NULL) {
    log_info("Fitting Phenotype model for DeepRVAT")
    all_res_list = list()
    all_models_list = list()
    for (gene_list in gene_lists) {
        log_info("Fitting model for gene_list: {gene_list}")
        all_models_list[[gene_list]] = list()
        estimate_list = list()
        # this_genes_to_keep = ifelse(is.null(genes_to_keep), NULL, list(genes_to_keep[[gene_list]]))

        if (is.null(genes_to_keep)) {
            this_genes_to_keep = NULL
        } else {
            this_genes_to_keep = list(genes_to_keep[[gene_list]])
        }
        log_info("computing avareaged burdens across n repeats: {n_deeprvat_repeats}")
        all_data_list = list()
        for (deeprvat_repeat in seq(0, n_deeprvat_repeats - 1)) {

            all_model_data = LoadDataForModel(y_data, covariate_cols = covariate_cols, phenotype_col = phenotype_col, top_q = top_q, this_out_dir = this_out_dir, method = "deeprvat", key_1 = gene_list, key_2 = deeprvat_repeat, genes_to_keep = this_genes_to_keep)
            all_data_list[paste("repeat", deeprvat_repeat, sep = "_")] = list(all_model_data)
        }

        for (key in c("train", "test")) {
            all_model_data[["data_list"]][[key]][["model_data"]] = AverageDeepRVATBurdens(all_data_list, key)
        }
        log_info("fitting model for avareaged burdens")
        data_list = all_model_data[["data_list"]]
        all_x_cols = all_model_data[["all_x_cols"]]
        model_output = FitLinearModel(data_list, all_x_cols, model_name = paste("avg_burden_deeprvat", gene_list, sep = "-"))

        all_res_list = append(all_res_list, list(model_output[["res"]]))

    }

    all_res_list = FlattenDfList(all_res_list)
    return(all_res_list)
}

PredictPhenoBaselineLinear <- function(y_data, covariate_cols, phenotype_col, this_out_dir, top_q, gene_lists, vtypes = c("plof", "missense"), genes_to_keep = NULL) {
    log_info("Fitting Phenotype model for Baseline")
    all_res_list = list()
    all_models_list = list()
    for (vtype in vtypes) {
        all_models_list[[vtype]] = list()
        inner_list = list()
        for (gene_list in gene_lists) {
            # this_genes_to_keep = ifelse(is.null(genes_to_keep), NULL, list(genes_to_keep[[gene_list]]))
            if (is.null(genes_to_keep)) {
                this_genes_to_keep = NULL
            } else {
                this_genes_to_keep = list(genes_to_keep[[gene_list]])
            }
            all_model_data = LoadDataForModel(y_data, covariate_cols = covariate_cols, phenotype_col = phenotype_col, top_q = top_q, this_out_dir = this_out_dir, method = "baseline", key_1 = vtype, key_2 = gene_list, genes_to_keep = this_genes_to_keep)
            data_list = all_model_data[["data_list"]]
            all_x_cols = all_model_data[["all_x_cols"]]
            model_output = FitLinearModel(data_list, all_x_cols, model_name = paste("baseline", vtype, gene_list, sep = "-"))

            all_res_list = append(all_res_list, list(model_output[["res"]]))
            all_models_list[[vtype]][[gene_list]] = model_output[["model"]]

        }
        # names(inner_list) = gene_lists all_res_list = append(all_res_list, list(inner_list))
    }
    # model_name_save = paste('baseline', phenotype_col, top_q, rank, sep = '-') saveRDS(all_models_list, paste0(data_dir, '/markdown_data/logistic_models/', model_name_save, '.Rds'))
    all_res_list = FlattenDfList(all_res_list)
    return(all_res_list)
}

PredictPhenoCovariatesLinear <- function(y_data, covariate_cols, phenotype_col, this_out_dir, top_q, use_top_q = TRUE) {
    log_info("Fitting Phenotype model for covariates only")

    all_model_data = LoadDataForModel(y_data, covariate_cols = covariate_cols, phenotype_col = phenotype_col, top_q = top_q, this_out_dir = this_out_dir, method = "deeprvat", key_1 = "baseline_only", key_2 = "1", use_top_q = use_top_q)  #can basically use any combination of method, key_1, key_2 since the covariates are all the same

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
