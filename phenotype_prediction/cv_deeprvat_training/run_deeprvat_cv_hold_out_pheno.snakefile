
from snakemake.utils import Paramspace
from snakemake.utils import min_version
from genopheno.name_mappings import BTYPES_DICT, PLOF_CONSEQUENCES
import os

min_version("6.0")

from pathlib import Path
import pandas as pd

configfile: 'config.yaml'

debug_flag = config.get('debug', False)
debug = '--debug ' if debug_flag else ''


conda_check = 'conda info | grep "active environment"'
cuda_visible_devices = 'echo CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES'

DEEPRVAT_ANALYSIS_DIR = os.environ['DEEPRVAT_ANALYSIS_DIR']
DEEPRVAT_DIR = os.environ['DEEPRVAT_DIR']

pipeline_dir = f'{DEEPRVAT_ANALYSIS_DIR}/phenotype_predction/cv_deeprvat_training/'
py_pipeline = f'python {pipeline_dir}'
py_deeprvat = f'python {DEEPRVAT_DIR}/deeprvat'


wildcard_constraints:
    repeat="\d+"

cv_splits = 5
repeats_to_compare = [6]
phenotypes = config['phenotypes']
phenotypes = phenotypes.keys() if type(phenotypes) == dict else phenotypes
phenotypes_testing = phenotypes
#phenotypes_testing = ['Triglycerides']


association_testing_maf = config.get('association_testing_maf', 0.001)

n_bags = config['deeprvat']['training'].get('n_bags')
n_repeats = config['deeprvat'].get('n_repeats', 6)
n_burden_chunks = config['deeprvat'].get('n_burden_chunks', 1) if not debug_flag else 2
n_trials = config['deeprvat']['hyperparameter_optimization']['n_trials']
n_regression_chunks = config['deeprvat'].get('n_regression_chunks', 40) if not debug_flag else 2

config['seed_genes']['phenotypes'] = config['phenotypes']
config['deeprvat']['phenotypes'] = config['phenotypes']


##### TO BE CHANGED BY USER ######
data_dir = 'X' #the directory where #train_multipheno_cv has been run.
#We use the seed gene discoveries and training and association data sets prepared for each phenotype from there to
#avoid rerunning
# The difference compared to train_multipheno_cv.snakefile is that we train models leaving out one {phenotype} (the directory is named after the held-out phenotype)
# and then compute the burdens for this held-out phenotype using the trained model. 



rule all:
    input:
        expand("{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/eval/significant.parquet",
               cv_split = range(cv_splits), phenotype=phenotypes_testing, type = ['deeprvat']),
        expand("{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/eval/all_results.parquet",
               cv_split = range(cv_splits), phenotype=phenotypes_testing, type = ['deeprvat']),
        expand('{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/burdens_test/burdens.zarr',
                phenotype = phenotypes_testing, cv_split=range(cv_splits), type = ['deeprvat']),
        expand('{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/burdens/burdens.zarr',
                phenotype = phenotypes_testing, cv_split=range(cv_splits), type = ['deeprvat'])

module deeprvat_workflow:
    snakefile: 
        f"{DEEPRVAT_ANALYSIS_DIR}/phenotype_prediction/cv_deeprvat_training/training_association_testing_with_prefix.snakefile"
    prefix:
        'cv_split{cv_split}/deeprvat'
    config:
        config['deeprvat']




use rule * from deeprvat_workflow exclude config, choose_training_genes, evaluate, train_bagging, regress, compute_burdens, compute_burdens_test, best_bagging_run, cleanup_burden_cache, link_burdens, link_burdens_test, all_burdens, combine_regression_chunks  as deeprvat_*

rule evaluate:
    input:
        associations = expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/repeat_{repeat}/results/burden_associations.parquet',
                              repeat=range(n_repeats)),
        config = f"{data_dir}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/hpopt_config.yaml",
    output:
        "{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/eval/significant.parquet",
        "{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/eval/all_results.parquet"
    threads: 1
    shell:
        'deeprvat_evaluate '
        # + debug +
        '--use-seed-genes '
        '--correction-method FDR '
        '{input.associations} '
        '{input.config} '
        '{wildcards.phenotype}/deeprvat/eval'

rule all_regression:
    input:
        expand('{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations.parquet',
               phenotype=phenotypes_testing, type=['deeprvat'], repeat=range(n_repeats), cv_split = range(cv_splits)),

use rule combine_regression_chunks from deeprvat_workflow as deeprvat_combine_regression_chunks with:
    input:
        testing = expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/repeat_{{repeat}}/results/burden_associations_{chunk}.parquet', chunk=range(n_regression_chunks)),
        replication = expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/repeat_{{repeat}}/results/burden_associations_replication_{chunk}.parquet', chunk=range(n_regression_chunks))
    output:
        testing = '{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations.parquet',
        replication = '{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations_replication.parquet'

use rule all_training_dataset from deeprvat_workflow as deeprvat_all_training_dataset with:
    input:
        input_tensor = expand(f'{data_dir}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/input_tensor.zarr',
                              phenotype=phenotypes,
                              cv_split=range(cv_splits)),
        covariates = expand(f'{data_dir}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/covariates.zarr',
                            phenotype=phenotypes,
                            cv_split=range(cv_splits)),
        y = expand(f'{data_dir}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/y.zarr',
                   phenotype=phenotypes,
                   cv_split=range(cv_splits)),


use rule regress from deeprvat_workflow as deeprvat_regress with:
    input:
        config = f"{data_dir}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/hpopt_config.yaml",
        chunks = expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/burdens/chunk{chunk}.finished',
                    chunk = range(n_burden_chunks))
    output:
        '{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations_{chunk}.parquet',
        '{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations_replication_{chunk}.parquet'
    params:
        prefix = '{phenotype}/cv_split{cv_split}/deeprvat'


rule all_burdens_zarr:
    input:
        expand('{p}/cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens_test/chunk{c}.finished',
                p = phenotypes_testing, c = range(n_burden_chunks), cv_split=range(cv_splits)), 
        expand('{p}/cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens/chunk{c}.finished',
                p =phenotypes_testing, c = range(n_burden_chunks), cv_split=range(cv_splits)),
    output:
        expand('{p}/cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens_test/burdens.zarr',
                p = phenotypes_testing, cv_split=range(cv_splits)),
        expand('{p}/cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens/burdens.zarr',
                p = phenotypes_testing, cv_split=range(cv_splits)),



use rule compute_burdens from deeprvat_workflow as deeprvat_compute_burdens with:
    input:
        checkpoints = expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/models/repeat_{repeat}/best/bag_{bag}.ckpt',
                             bag=range(n_bags), repeat=range(n_repeats)),
        dataset = f'{data_dir}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/association_dataset.pkl',
        config = '{phenotype}/cv_split{cv_split}/deeprvat/models/repeat_0/trial0/config.yaml', #TODO make this more generic #TODO check if this is correct
        reversed = "{phenotype}/cv_split{cv_split}/deeprvat/models/reverse_finished.tmp",
    output:
        '{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/burdens/chunk{chunk}.finished'
    params:
        prefix= '{phenotype}/cv_split{cv_split}/deeprvat'

use rule reverse_models from deeprvat_workflow as deeprvat_reverse_models with:
    input:
        checkpoints = expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/models/repeat_{repeat}/best/bag_{bag}.ckpt',
                             bag=range(n_bags), repeat=range(n_repeats)),
        config = '{phenotype}/cv_split{cv_split}/deeprvat/models/repeat_0/config.yaml', 
    output:
        "{phenotype}/cv_split{cv_split}/deeprvat/models/reverse_finished.tmp"

use rule best_training_run from deeprvat_workflow as deeprvat_best_training_run with:
    input:
        expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/models/repeat_{{repeat}}/trial{trial_number}/config.yaml',
               trial_number=range(n_trials)),
    output:
        checkpoints = expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/models/repeat_{{repeat}}/best/bag_{bag}.ckpt',
                             bag=range(n_bags)),
        config = '{phenotype}/cv_split{cv_split}/deeprvat/models/repeat_{repeat}/config.yaml'
    params:
        prefix = '{phenotype}/cv_split{cv_split}/deeprvat'

rule all_training:
    input:
        expand('{phenotype}/cv_split{cv_split}/deeprvat/models/repeat_{repeat}/trial{trial_number}/config.yaml',
               phenotype = phenotypes, 
               cv_split = range(cv_splits),
               repeat=range(n_repeats),
               trial_number=range(n_trials)),
        expand('{phenotype}/cv_split{cv_split}/deeprvat/models/repeat_{repeat}/trial{trial_number}/finished.tmp',
               phenotype = phenotypes, 
               cv_split = range(cv_splits),
               repeat=range(n_repeats),
               trial_number=range(n_trials))


use rule train from deeprvat_workflow as deeprvat_train with:
    input:
        config = f'{data_dir}/cv_split{{cv_split}}/deeprvat/config.yaml',
        input_tensor = expand(f'{data_dir}/cv_split{{{{cv_split}}}}/deeprvat/{{phenotype}}/deeprvat/input_tensor.zarr',
                              phenotype=phenotypes),
        covariates = expand(f'{data_dir}/cv_split{{{{cv_split}}}}/deeprvat/{{phenotype}}/deeprvat/covariates.zarr',
                            phenotype=phenotypes),
        y = expand(f'{data_dir}/cv_split{{{{cv_split}}}}/deeprvat/{{phenotype}}/deeprvat/y.zarr',
                   phenotype=phenotypes),
    output:
        config = '{phenotype}/cv_split{cv_split}/deeprvat/models/repeat_{repeat}/trial{trial_number}/config.yaml',
        finished = '{phenotype}/cv_split{cv_split}/deeprvat/models/repeat_{repeat}/trial{trial_number}/finished.tmp'
    resources:
        mem_mb = 2000000,        # Using this value will tell our modified lsf.profile not to set a memory resource
        load = 8000,
        gpus = 1
    params:
        prefix = '{phenotype}/cv_split{cv_split}/deeprvat',
        phenotypes = lambda wildcards: " ".join( #TODO like need the prefix here as well
            [f"--phenotype {p} "
             f"{data_dir}/cv_split{wildcards.cv_split}/deeprvat/{p}/deeprvat/input_tensor.zarr "
             f"{data_dir}/cv_split{wildcards.cv_split}/deeprvat/{p}/deeprvat/covariates.zarr "
             f"{data_dir}/cv_split{wildcards.cv_split}/deeprvat/{p}/deeprvat/y.zarr"
             for p in phenotypes if p != wildcards.phenotype]) 



############################### Computation of test set deeprvat burdens #####################################################
############################################################################################################################


use rule compute_burdens_test from deeprvat_workflow as deeprvat_compute_burdens_test with:
    input:
        checkpoints = expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/models/repeat_{repeat}/best/bag_{bag}.ckpt',
                             bag=range(n_bags), repeat=range(n_repeats)),
        dataset = f'{data_dir}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/association_dataset_test.pkl',
        config = '{phenotype}/cv_split{cv_split}/deeprvat/models/repeat_0/trial0/config.yaml', 
    output:
        '{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/burdens_test/chunk{chunk}.finished'
    params:
        prefix = '{phenotype}/cv_split{cv_split}/deeprvat',


rule compute_burdens_test:
    priority: 10
    input:
        checkpoints = expand('{{phenotype}}/cv_split{{cv_split}}/deeprvat/models/repeat_{repeat}/best/bag_{bag}.ckpt',
                             bag=range(n_bags), repeat=range(n_repeats)),
        dataset = f'{data_dir}/cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/association_dataset_test.pkl',
        config = '{phenotype}/cv_split{cv_split}/deeprvat/models/repeat_0/trial0/config.yaml', 
    output:
        '{phenotype}/cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/burdens_test/chunk{chunk}.finished'
    threads: 8
    shell:
        ' && '.join([
            ('deeprvat_associate compute-burdens '
             + debug +
             ' --n-chunks '+ str(n_burden_chunks) + ' '
             '--chunk {wildcards.chunk} '
             '--dataset-file {input.dataset} '
             '{input.config} '
             '{input.checkpoints} '
             '{wildcards.phenotype}/deeprvat/burdens_test'),
            'touch {output}'
        ])

