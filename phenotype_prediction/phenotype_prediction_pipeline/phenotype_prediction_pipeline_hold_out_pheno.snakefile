from snakemake.utils import Paramspace
from snakemake.utils import min_version
import os
import copy
min_version("6.0")

configfile: 'config_eval.yaml'
from pathlib import Path
import pandas as pd

# configfile: 'config.yaml'
conda_check = 'conda info | grep "active environment"'


cv_splits = config.get('cv_splits', 5)



phenotypes = config.get('phenotypes')
# phenotypes = phenotypes[:2]
phenotypes = [phenotype.replace('_standardized', '') for phenotype in phenotypes]
#phenotypes = ['Triglycerides']

DEEPRVAT_ANALYSIS_DIR = os.environ['DEEPRVAT_ANALYSIS_DIR']
code_dir = f'{DEEPRVAT_ANALYSIS_DIR}/phenotype_prediction/phenotype_prediction_pipeline/' 

top_bottom_quantiles = ['topq', 'bottomq']

regression_config = config['regression']
phenotype_suffixes = regression_config['pheno_suffixes']
top_quantiles = regression_config['top_quantiles']
# fdrs = regression_config['fdrs']
fdrs = regression_config.get('fdrs', [0.05])
simple_btypes = copy.deepcopy(regression_config['r_config']['vtypes'])
btypes = copy.deepcopy(simple_btypes)
btypes.append('deeprvat')
burden_btypes = [btype for btype in btypes if '_' not in btype]

##### TO BE CHANGED BY USER ######
data_dir = 'X' #the directory where #train_multipheno_cv has been run.

py = f'python {code_dir}/'



rule all_prs_prediction_hold_out:
    input:
        expand('phenotype_prediction/models/linear_models/plotting_data/all_recomputed_metrics_test_{phenotype_suffix}_{fdr}.Rds',
                phenotype_suffix = phenotype_suffixes,
                fdr = fdrs),
        expand('phenotype_prediction/models/linear_models/plotting_data/plot_df_list_{phenotype_suffix}_{fdr}.Rds',
                phenotype_suffix = phenotype_suffixes,
                fdr = fdrs),
        expand('phenotype_prediction/models/logistic_models/plotting_data/{logistic_data}_{phenotype_suffix}_{fdr}.Rds',
                phenotype_suffix = phenotype_suffixes,
                fdr = fdrs,
                logistic_data = [#'combined_res',
                 'combined_metrics'])

module prs_pipeline:
    snakefile: 
        "phenotype_prediction_pipeline.snakefile"
    prefix:
        './'
    config:
        config

use rule * from prs_pipeline exclude prep_data_for_r, extract_burdens as prs_*



use rule make_plot_df_list_data_linear_model from prs_pipeline as prs_make_plot_df_list_data_linear_model with:
    input:
        expand('phenotype_prediction/models/linear_models/{phenotype}/{{phenotype_suffix}}_fdr-{{fdr}}.Rds', 
                phenotype = phenotypes)
    output: 'phenotype_prediction/models/linear_models/plotting_data/plot_df_list_{phenotype_suffix}_{fdr}.Rds'

use rule prep_metrics_linear_model from prs_pipeline as prs_prep_metrics_linear_model with:
    input:
        expand('phenotype_prediction/models/linear_models/{phenotype}/{{phenotype_suffix}}_fdr-{{fdr}}.Rds', 
                phenotype = phenotypes)
    output: 'phenotype_prediction/models/linear_models/plotting_data/all_recomputed_metrics_test_{phenotype_suffix}_{fdr}.Rds'

use rule prep_metrics_logistic_model from prs_pipeline as prs_prep_metrics_logistic_model with:
    input:
        expand('phenotype_prediction/models/logistic_models/{phenotype}/{{phenotype_suffix}}_{top_bottom_q}-{top_q}_fdr-{{fdr}}.Rds', 
                phenotype = phenotypes,
                top_q = top_quantiles,
                top_bottom_q = top_bottom_quantiles)
    output: 
        metrics = 'phenotype_prediction/models/logistic_models/plotting_data/combined_metrics_{phenotype_suffix}_{fdr}.Rds',

use rule fit_logistic_regression_model_fdr from  prs_pipeline as prs_fit_logistic_regression_model_fdr with:
    input: 
        expand('phenotype_prediction/r_data/{{phenotype}}/cv_split{cvsplit}/data.finished',
                    cvsplit = range(cv_splits))
    output: 
        'phenotype_prediction/models/logistic_models/{phenotype}/{phenotype_suffix}_{top_bottom_q}-{top_q}_fdr-{fdr}.Rds'


use rule fit_linear_regression_model_fdr from  prs_pipeline as prs_fit_linear_regression_model_fdr with:
    input: 
        expand('phenotype_prediction/r_data/{{phenotype}}/cv_split{cvsplit}/data.finished',
                    cvsplit = range(cv_splits))
    output: 
        'phenotype_prediction/models/linear_models/{phenotype}/{phenotype_suffix}_fdr-{fdr}.Rds'


use rule prep_data_for_r from prs_pipeline as prs_prep_data_for_r with:
    input:
        burdens = expand('phenotype_prediction/burdens/{{phenotype}}/cv_split{cvsplit}/{split}_{btype}_burdens.pickle', 
                split = ['train', 'test'],
                cvsplit = range(cv_splits),
                btype = btypes),
        config = 'config_eval.yaml',
        dataset_train = f'{data_dir}/cv_split{{cvsplit}}/deeprvat/{{phenotype}}/deeprvat/association_dataset.pkl',
        dataset_test = f'{data_dir}/cv_split{{cvsplit}}/deeprvat/{{phenotype}}/deeprvat/association_dataset_test.pkl',
        discoveries = '{phenotype}/cv_split{cvsplit}//deeprvat/{phenotype}/deeprvat/eval/significant.parquet', 
    output:
        'phenotype_prediction/r_data/{phenotype}/cv_split{cvsplit}/data.finished'
    params:
        dataset_dir = f'{data_dir}/cv_split{{cvsplit}}/deeprvat/{{phenotype}}/deeprvat/',
        burdens_dir = 'phenotype_prediction/burdens/{phenotype}/cv_split{cvsplit}',
        out_dir = 'phenotype_prediction/r_data/{phenotype}/cv_split{cvsplit}',
        burden_types = [f'--burden-type {b} 'for b in btypes]


# # for simple burdens use multipheno_cv experiments because they stay the same


burden_dir_mapper = {
    btype: "baseline_burdens" if btype != "deeprvat" else btype
    for btype in burden_btypes
}
split_mapper = {'train': '', 'test' :'_test'}

use rule extract_burdens from prs_pipeline as prs_extract_burdens with:
    input:
        #discoveries don't have to be combined since we anyways only run association testing for the held-out phenotyp (so there is only one file)
        discoveries = '{phenotype}/cv_split{cvsplit}//deeprvat/{phenotype}/deeprvat/eval/significant.parquet', 
        burdens = lambda wildcards: f'{wildcards.phenotype}/cv_split{wildcards.cvsplit}/{burden_dir_mapper[wildcards.btype]}/{wildcards.phenotype}/{wildcards.btype}/burdens{split_mapper[wildcards.split]}/burdens.zarr' \
            if wildcards.btype == 'deeprvat' else f'{data_dir}/cv_split{wildcards.cvsplit}/{burden_dir_mapper[wildcards.btype]}/{wildcards.phenotype}/{wildcards.btype}/burdens{split_mapper[wildcards.split]}/burdens.zarr' 
    output:
        'phenotype_prediction/burdens/{phenotype}/cv_split{cvsplit}/{split}_{btype}_burdens.pickle'
    resources:
        mem_mb = 16000,
        load = 16000
    params:
        burden_input_dir = lambda wildcards: f'{wildcards.phenotype}/cv_split{wildcards.cvsplit}/{burden_dir_mapper[wildcards.btype]}/{wildcards.phenotype}/{wildcards.btype}/burdens{split_mapper[wildcards.split]}' \
            if wildcards.btype == 'deeprvat' else f'{data_dir}/cv_split{wildcards.cvsplit}/{burden_dir_mapper[wildcards.btype]}/{wildcards.phenotype}/{wildcards.btype}/burdens{split_mapper[wildcards.split]}' 



include: '../cv_deeprvat_training/run_deeprvat_cv_hold_out_pheno.snakefile' 
#when including this, comment out # include: '../cv_deeprvat_training/run_deeprvat_cv.snakefile' in run_deeprvat_cv_hold_out_pheno.snakefile
