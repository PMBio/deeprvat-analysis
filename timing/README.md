# Get compute times for different RVAT methods

In this directory, we run the experiments and evaluation required to get Supp. Fig. 1.1 and Figure 3c. 

For this, association testing is run with all methods (Burden/SKAT, Monti, STAAR, DeepRVAT) in the experiment directories provided here using the relevant snakfiles (see `comparison_methos` and `association_testing`). 
The association testing is only run on a random set of 1000 genes and compute times are averaged across them. 


For DeepRVAT, we only run the association testing here (`pretrained_models` have to be linked from `association_testing/paper_experiment`) using the [DeepRVAT association testing pipeline](https://github.com/PMBio/deeprvat/blob/main/pipelines/association_testing_pretrained.snakefile). 
The DeepRVAT training and burden computation times have been retrieved from the log files after running the `association_testing/paper_experiment`. 





