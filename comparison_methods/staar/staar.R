library(STAAR)
library(rhdf5)
library(dplyr)
library(arrow)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)


genotype_file <- args[1]
annotation_file <- args[2]
x_file <- args[3]
y_file <- args[4]
out_path = args[5]
min_eac = 10  #as in original STAAR paper

print(paste("Using min EAC for variants per gene of", min_eac, sep = " "))
print(out_path)

covariates = h5read(x_file, "X")
covariates_dt = as.data.frame(t(covariates))

# Constant column not expected by null model fitting causes error when running STAAR and therefore has to be removed
if (length(unique(covariates_dt[, 1])) == 1) {
    print("removing constant column from covariates data frame")
    covariates_dt = covariates_dt[, seq(2, ncol(covariates_dt))]
}

y_data = h5read(y_file, "y")
y_dt = as.data.frame(t(y_data))
colnames(y_dt) = c("Y")
print(length(unique(y_dt)))
print(length(unique(y_dt)) == 2)
if (length(unique(y_dt)) == 2) {
    print("Binary y, using binomial family for glm")
    glm_family = "binomial"
} else {
    print("Continuous y, using binomial family for glm")
    glm_family = "gaussian"
}

null_model_data = cbind(y_dt, covariates_dt)

null_model_formula = paste("Y ~ ", paste(colnames(covariates_dt), collapse = "+"))
print(null_model_formula)

obj_nullmodel <- fit_null_glm(null_model_formula, data = null_model_data, family = glm_family)

genes <- h5ls(genotype_file)$name
print(paste("Running STAAR on ", length(genes), " genes"))

if (!all(sort(h5ls(annotation_file)$name) == sort(genes))) {
    stop("Gene lists in HDF5 files don't agree")
}

runStaarOnGene = function(this_gene, genotype_file, annotation_file, obj_nullmodel) {
    this_geno = t(h5read(genotype_file, this_gene))
    EAC = sum(this_geno)
    this_anno = t(h5read(annotation_file, this_gene))
    print(dim(this_geno))
    print(dim(this_anno))
    print(paste("EAC", EAC))

    if ((nrow(this_anno) > 10000) | (EAC < min_eac)) {
        print("Not running STAAR because of > 10000 variants or EAC < min_eac")
        print("Setting pval to NA")
        return(c(pval = NA, EAC = EAC, time = NA))
    } else {

        tryCatch({
            start_time = Sys.time()
            staar_res = STAAR(genotype = this_geno, obj_nullmodel = obj_nullmodel, annotation_phred = this_anno, rare_maf_cutoff = 0.001)
            end_time = Sys.time()
            pval = staar_res$results_STAAR_O
            time_diff = end_time - start_time
            time_diff = str_extract(time_diff, "\\d+\\.*\\d*")
            return(c(pval = pval, EAC = EAC, time = time_diff))
        }, error = function(error_message) {
            message("STAAR error occured")
            message("And below is the error message from R:")
            message(error_message)
            return(c(pval = NA, EAC = EAC, time = NA))
        })
    }
}


pval_list = list()
eac_list = list()
time_list_outer = list()
time_list_inner = list()
for (i in seq(length(genes))) {

    this_gene = genes[[i]]
    print(paste(Sys.time(), "Running STAAR for gene ", i, this_gene, " out of ", length(genes), "genes"))
    staar_res = runStaarOnGene(this_gene, genotype_file, annotation_file, obj_nullmodel)
    pval_list[[this_gene]] = staar_res[["pval"]]
    eac_list[[this_gene]] = staar_res[["EAC"]]
    time_list_inner[[this_gene]] = staar_res[["time"]]


}
res = tibble(gene = names(pval_list), pval = unlist(pval_list), EAC = unlist(eac_list), time_inner = unlist(time_list_inner))
print(head(res))
saveRDS(res, out_path)
print("finished")
