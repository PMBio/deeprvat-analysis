# used in monit_replication.R and staar_replication.R
aggPvalsToGene = function(df, combine_pvals, grouping_cols) {
    if (tolower(combine_pvals) == "cct") {
        print("Combining pvals for gene-trait combi using CCT")
        df = df %>%
            drop_na(pval) %>%
            group_by(across(all_of(grouping_cols))) %>%
            summarize(pval = CCT(pval)) %>%
            ungroup() %>%
            mutate(pval_agg_method = "cct_combined")

    } else if (tolower(combine_pvals) == "bonferroni") {
        print("Combining pvals for gene-trait combi using bonferroni correction and min pval")
        df = df %>%
            group_by(across(all_of(grouping_cols))) %>%
            mutate(pval = p.adjust(pval, method = "bonferroni")) %>%
            summarize(pval = min(pval)) %>%
            ungroup() %>%
            mutate(pval_agg_method = "bonferroni_min_combined")
    } else {
        warning("invalid pvalue combination method")
    }
    return(df %>%
        ungroup())
}

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
