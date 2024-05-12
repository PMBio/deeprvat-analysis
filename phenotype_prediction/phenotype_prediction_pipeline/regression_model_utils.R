library(tibble)
library(arrow)
library(groupdata2)
library(logger)
### name and input path settings ##########


genetic_pcs = paste0("genetic_PC_", 1:20)
prs = "prs"
age = "age"
sex = "genetic_sex"

covariates_dict <- list(c(prs, age, genetic_pcs), c(prs, genetic_pcs), c(genetic_pcs), c(prs, age, sex, genetic_pcs), c(prs, genetic_pcs), c(genetic_pcs), c(prs, age, sex, genetic_pcs), c(prs, genetic_pcs), c(genetic_pcs), c(prs, age, sex, genetic_pcs), c(prs, genetic_pcs), c(genetic_pcs))


names(covariates_dict) = c(paste0("standardized"), paste0("standardized_resid_wo_prs"), paste0("standardized_resid_with_prs"), paste0("orig"), paste0("resid_wo_prs"), paste0("resid_with_prs"), paste0("irnt"), paste0("irnt_resid_wo_prs"), paste0("irnt_resid_with_prs"), paste0("qt"),
    paste0("qt_resid_wo_prs"), paste0("qt_resid_with_prs"))

# TODO delete this either here or in regression_analyis.R
model_names_ordered = c("covariates", "covariates-PRS", "baseline-missense-baseline_only", "baseline-missense-deeprvat_discoveries", "baseline-missense-genebass", "baseline-missense-deeprvat_all_samples", "baseline-plof-baseline_only", "baseline-plof-deeprvat_discoveries", "baseline-plof-genebass",
    "baseline-plof-deeprvat_all_samples", "deeprvat-baseline_only", "deeprvat-deeprvat_discoveries", "deeprvat-genebass", "deeprvat-deeprvat_all_samples")


model_names_ordered_renamed = c("covariates", "covariates-PRS", "missense-baseline_discoveries", "missense-deeprvat_discoveries", "missense-genebass", "missense-deeprvat_discoveries_all_samples", "plof-baseline_discoveries", "plof-deeprvat_discoveries", "plof-genebass", "plof-deeprvat_discoveries_all_samples",
    "deeprvat-baseline_discoveries", "deeprvat-deeprvat_discoveries", "deeprvat-genebass", "deeprvat-deeprvat_discoveries_all_samples")

gene_list_names = c(NA, NA, "baseline_discoveries", "deeprvat_discoveries", "genebass", "deeprvat_discoveries_all_samples", "baseline_discoveries", "deeprvat_discoveries", "genebass", "deeprvat_discoveries_all_samples", "baseline_discoveries", "deeprvat_discoveries", "genebass", "deeprvat_discoveries_all_samples")

##########################################



GetColumnNames <- function(phenotype_col_suffix, phenotype, covariates_dict) {

    if (phenotype_col_suffix == "orig") {
        phenotype_col = phenotype
        covariate_cols = unlist(covariates_dict["orig"])
    } else {
        phenotype_col = paste0(phenotype, "_", phenotype_col_suffix)
        covariate_cols = unlist(covariates_dict[phenotype_col_suffix])
    }
    phenotype_col

    return(list(covariate_cols = covariate_cols, phenotype_col = phenotype_col))
}


LoadYData <- function(this_out_dir, splits = c("train", "test")) {
    y_data = list()
    for (split in splits) {
        this_y_data = read_parquet(paste0(this_out_dir, "/", split, "-y.parquet")) %>%
            column_to_rownames("sample")
        y_data = append(y_data, list(this_y_data))
    }
    names(y_data) = splits
    return(y_data)
}

LoadDataForModel <- function(y_data, covariate_cols, phenotype_col, top_q, this_out_dir, btype, gene_list, splits = c("train", "test"), genes_to_keep = NULL, use_top_q = TRUE) {


    data_list = list()
    for (split in splits) {
        x_file_path <- paste0(this_out_dir, "/", split, "-", btype, "-", gene_list, "-", "x.parquet")
        print(x_file_path)
        this_x_data = read_parquet(x_file_path)
        if (split == "train") {
            gene_cols <- colnames(this_x_data)[grep("^gene_", colnames(this_x_data))]
            if (is.null(genes_to_keep)) {
                log_info("Using all genes in input data frame")
                # gene_cols <- colnames(this_x_data)[grep('^gene_', colnames(this_x_data))]
                all_x_cols = c(covariate_cols, gene_cols)
            } else {
                log_info("Only using genes in genes_to_keep")
                if (length(unlist(genes_to_keep)) > 0) {
                  # ifelse allows to use now genes (i.e. the covariates only model)
                  gene_cols_select = paste("gene", unlist(genes_to_keep), sep = "_")
                  gene_cols = Reduce(intersect, list(gene_cols, gene_cols_select))

                  all_x_cols = c(covariate_cols, gene_cols)
                } else {
                  all_x_cols = covariate_cols
                }
            }
            # all_x_cols = c(covariate_cols, gene_cols)
        }
        X = this_x_data[all_x_cols]
        # print(colnames(y_data[[split]])) print(phenotype_col)
        Y = y_data[[split]][phenotype_col]
        if (use_top_q) {
            log_info("Using top quantile")
            thres = quantile(Y[[phenotype_col]], 1 - top_q)
            Y_bin = Y > thres
        } else {
            log_info("Using bottom quantile")
            thres = quantile(Y[[phenotype_col]], top_q)
            Y_bin = Y < thres
        }
        # assert all(as.integer(rownames(Y_bin)) == this_x_data$sample)
        model_data = cbind(X, as_tibble_col(Y_bin, column_name = "y_bin"))


        this_data_list = list(X = X, Y = Y, Y_bin = Y_bin, model_data = model_data)

        data_list = append(data_list, list(this_data_list))
    }
    names(data_list) = splits
    return(list(data_list = data_list, all_x_cols = all_x_cols))
}

EvalLogisticModel <- function(fitted_results, truth, y_true_cont, plot_title = "", decision_threshold = 0.5) {
    fitted_results_bin <- ifelse(fitted_results > decision_threshold, 1, 0)
    this_res = tibble(estimate = fitted_results, estimate_bin = as.factor(fitted_results_bin), truth = as.factor(truth), Y = y_true_cont)

    class_and_probs_metrics <- metric_set(roc_auc, pr_auc)

    this_metrics = this_res %>%
        class_and_probs_metrics(truth, estimate, estimate = estimate_bin, event_level = "second")
    auprc = this_metrics %>%
        filter(str_detect(.metric, "pr_auc")) %>%
        select(.estimate)
    # this_plot = autoplot(pr_curve(this_res, truth, estimate, event_level = 'second')) + ggtitle(paste0(plot_title, ', AUPRC: ', round(auprc, 4)*100, '%')) show(this_plot)

    return(list(res = this_res, metrics = this_metrics))
}

FitLogisticModel <- function(data_list, all_x_cols, upsample = TRUE, use_weight = FALSE, model_name = "", decision_threshold = 0.5) {
    log_info("Fitting model")
    model_formula = paste0("y_bin ~ ", paste(all_x_cols, collapse = " + "))
    # print(paste0('Model formula: ', model_formula))
    model_data = data_list[["train"]][["model_data"]]
    if (upsample) {
        model_data = upsample(model_data, "y_bin")
        # print(paste0('Upsampling data ', nrow(model_data)))

    }
    model = glm(model_formula, data = model_data, family = binomial(link = "logit"), weights = NULL)

    test_split = "test"
    fitted_results = predict(model, newdata = data_list[[test_split]][["model_data"]], type = "response")

    print(summary(model))

    combined_res = EvalLogisticModel(fitted_results, as.integer(data_list[[test_split]][["Y_bin"]]), data_list[[test_split]][["Y"]] %>%
        pull(var = 1), model_name)

    combined_res = mapply(cbind, combined_res, model_name = model_name, SIMPLIFY = F)
    return_data = list(res = combined_res, model = model$coefficients)
    # return(combined_res)
    return(return_data)

}


PredictPhenoRareBurden <- function(y_data, covariate_cols, phenotype_col, this_out_dir, top_q, gene_lists, btypes = c("plof", "missense"), genes_to_keep = NULL, use_top_q = TRUE) {
    log_info("Fitting Phenotype model for Rare burden")
    all_res_list = list()
    all_models_list = list()
    for (btype in btypes) {
        all_models_list[[btype]] = list()
        inner_list = list()
        for (gene_list in gene_lists) {
            # this_genes_to_keep = ifelse(is.null(genes_to_keep), NULL, list(genes_to_keep[[gene_list]]))
            if (is.null(genes_to_keep)) {
                this_genes_to_keep = NULL
            } else {
                this_genes_to_keep = list(genes_to_keep[[gene_list]])
            }
            all_model_data = LoadDataForModel(y_data, covariate_cols = covariate_cols, phenotype_col = phenotype_col, top_q = top_q, this_out_dir = this_out_dir, btype = btype, gene_list = gene_list, genes_to_keep = this_genes_to_keep, use_top_q = use_top_q)
            data_list = all_model_data[["data_list"]]
            all_x_cols = all_model_data[["all_x_cols"]]
            model_output = FitLogisticModel(data_list, all_x_cols, use_weight = FALSE, model_name = paste(btype, gene_list, sep = "-"))

            all_res_list = append(all_res_list, list(model_output[["res"]]))
            all_models_list[[btype]][[gene_list]] = model_output[["model"]]

        }
    }
    all_res_list = FlattenDfList(all_res_list)
    return(all_res_list)
}


PredictPhenoCovariates <- function(y_data, covariate_cols, phenotype_col, this_out_dir, top_q, use_top_q = TRUE) {
    log_info("Fitting Phenotype model for covariates only")
    # use any btype here just to get the covariates
    all_model_data = LoadDataForModel(y_data, covariate_cols = covariate_cols, phenotype_col = phenotype_col, top_q = top_q, this_out_dir = this_out_dir, btype = "deeprvat", gene_list = "baseline_only", use_top_q = use_top_q)  #can basically use any combination of method, key_1, key_2 since the covariates are all the same

    data_list = all_model_data[["data_list"]]
    all_x_cols = all_model_data[["all_x_cols"]]
    covariates_to_keep = c(age, sex, genetic_pcs)
    cols_to_keep = unlist(to_list(for (string in all_x_cols) if (any(sapply(covariates_to_keep, grepl, string)))
        string))

    all_res_list = list()
    for (useprs in c(FALSE, TRUE)) {
        log_info("use prs {useprs}")
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
        model_output = FitLogisticModel(data_list, this_cols_to_keep, use_weight = FALSE, model_name = model_name)
        all_res_list = append(all_res_list, list(model_output[["res"]]))

        # names(inner_list) = gene_lists all_res_list = append(all_res_list, list(inner_list))
    }
    all_res_list = FlattenDfList(all_res_list)
    return(all_res_list)
}


library(comprehenr)
FlattenDfList <- function(df_list) {
    df_names = names(df_list[[1]])
    concat_res_list = list()
    for (name in df_names) {
        this_dfs = to_list(for (i in seq(length(df_list))) df_list[[i]][[name]])
        res_concat = do.call(rbind, this_dfs)
        concat_res_list = append(concat_res_list, list(res_concat))
    }
    names(concat_res_list) = df_names
    return(concat_res_list)
}


