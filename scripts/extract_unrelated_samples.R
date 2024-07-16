library("ukbtools")
library(dplyr)
library(ggplot2)
library(arrow)
library(stringr)
library(tidyr)
library(reticulate)



phenotype_df = read_parquet('phenotypes.parquet') #as in https://github.com/PMBio/deeprvat/tree/main/example
all_samples = phenotype_df %>% pull('sample')


ukb_relatedness = read.csv('ukb_rel_a81358_s488120.dat', sep = ' ') #Kinship file provided by UK Biobank (Resource 668)
ukb_relatedness %>% arrange(Kinship)


print(length(all_samples))


print('extracting unrelated samples')

samples_to_remove = ukb_gen_samples_to_remove(ukb_relatedness, ukb_with_data = all_samples, 
                                                        cutoff = 0.0884) # includes pairs with greater than 3rd-degree relatedness
print(length(samples_to_remove))
unrelated_samples_to_keep = setdiff(all_samples, samples_to_remove)


print('writing')

unrelated_samples_out_file = '3rd_degree_unrelated_samples.pkl'
py_save_object(unrelated_samples_to_keep, unrelated_samples_out_file, pickle = "pickle")

