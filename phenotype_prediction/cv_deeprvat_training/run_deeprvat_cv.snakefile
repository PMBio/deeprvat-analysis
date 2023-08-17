
from snakemake.utils import Paramspace
from snakemake.utils import min_version
from name_mappings import BTYPES_DICT, PLOF_CONSEQUENCES
import os
min_version("6.0")

from pathlib import Path
import pandas as pd

configfile: 'config.yaml'

debug_flag = config.get('debug', False)
# debug_flag = True #TODO change this
debug = '--debug ' if debug_flag else ''



conda_check = 'conda info | grep "active environment"'
cuda_visible_devices = 'echo CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES'

DEEPRVAT_ANALYSIS_DIR = os.environ['DEEPRVAT_ANALYSIS_DIR']
DEEPRVAT_DIR = os.environ['DEEPRVAT_DIR']

pipeline_dir = f'{DEEPRVAT_ANALYSIS_DIR}/phenotype_prediction/cv_deeprvat_training/'
py_pipeline = f'python {pipeline_dir}'
py_deeprvat = f'python {DEEPRVAT_DIR}/deeprvat'

wildcard_constraints:
    repeat="\d+"

cv_splits = 5
alt_burdens_chunks = 30
repeats_to_compare = [6]
phenotypes = config['phenotypes']
phenotypes = phenotypes.keys() if type(phenotypes) == dict else phenotypes
burden_phenotype = phenotypes_testing[0]


association_testing_maf = config.get('association_testing_maf', 0.001)

n_bags = config['deeprvat']['training'].get('n_bags')
n_repeats = config['deeprvat'].get('n_repeats', 6)
n_burden_chunks = config['deeprvat'].get('n_burden_chunks', 1) if not debug_flag else 2
n_regression_chunks = config['deeprvat'].get('n_regression_chunks', 40) if not debug_flag else 2
n_regression_chunks = 4

config['seed_genes']['phenotypes'] = config['phenotypes']
config['deeprvat']['phenotypes'] = config['phenotypes']


btypes = config['alternative_burdens']['alt_burdens_data']['dataset_config']['rare_embedding']['config']['annotations']
if 'Consequence_stop_lost' in btypes:
    btypes.append('is_plof')

btypes = list(set(btypes) - set(PLOF_CONSEQUENCES))

btypes = [BTYPES_DICT[btype] if btype in BTYPES_DICT else btype for btype in btypes]
btypes.append('plof')


rule all:
    input:
        expand("cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/eval/significant.parquet",
               cv_split = range(cv_splits), phenotype=phenotypes_testing),
        expand("cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/eval/all_results.parquet",
               cv_split = range(cv_splits), phenotype=phenotypes_testing),
        expand('cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens_test/burdens.zarr',
                p = phenotypes_testing[0], cv_split=range(cv_splits)),
        expand('cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens/burdens.zarr',
                p = phenotypes_testing[0], cv_split=range(cv_splits)),
        expand('cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens_test/chunk{c}.linked',
                p = [pheno for pheno in phenotypes_testing if pheno != burden_phenotype],
                 c = range(n_burden_chunks), cv_split=range(cv_splits)),
        expand('cv_split{cv_split}/baseline/{p}/eval/burden_associations.parquet',
                p = phenotypes_testing, cv_split=range(cv_splits)),
        expand('cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/burdens/burdens.zarr',
                cv_split=range(cv_splits),
                btype = btypes,
                phenotype = phenotypes_testing),      
        expand('cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/burdens_test/burdens.zarr',
                cv_split=range(cv_splits),
                btype = btypes,
                phenotype = phenotypes_testing), 
        expand('cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/burdens/linking.finished',
                cv_split=range(cv_splits),
                btype = btypes,
                phenotype = phenotypes_testing),


rule spread_config:
    input:
        config = 'config.yaml'
    output:
        baseline = 'cv_split{cv_split}/baseline/config.yaml',
        deeprvat = 'cv_split{cv_split}/deeprvat/config.yaml',
        alternative_burdens = 'cv_split{cv_split}/alternative_burdens/config.yaml'
    params:
        out_path = 'cv_split{cv_split}/'
    threads: 1
    resources:
        mem_mb = 1024,
        load = 1000
    shell:
        ' && '.join([
            conda_check,
            py_pipeline + 'cv_utils.py spread-config '
            '-m seed_genes -m deeprvat -m alternative_burdens '
            '--fold {wildcards.cv_split} '
            ' {input.config} {params.out_path}'
        ])



rule all_config:
    input: 
        expand('cv_split{cv_split}/alternative_burdens/config.yaml', cv_split = range(cv_splits))


##### run seed gene discovery in every cv split  ######################################################################################
#########################################################################################################
# comment out if the baseline has already run
module seed_gene_workflow:
    snakefile: 
        f"{DEEPRVAT_DIR}/pipelines/seed_gene_discovery.snakefile"
    prefix: 
        'cv_split{cv_split}/baseline'
    config: 
        config['seed_genes']


use rule * from seed_gene_workflow exclude config as seed_gene_*

use rule config from seed_gene_workflow as seed_gene_config with:
    input:
        config = 'cv_split{cv_split}/baseline/config.yaml', 
    params:
        rare_maf = association_testing_maf
        

##########################################

# ############################### Run DeepRVAT ##############################################################
# ###########################################################################################################
module deeprvat_workflow:
    snakefile: 
        f"{DEEPRVAT_ANALYSIS_DIR}/phenotype_prediction/cv_deeprvat_training/training_association_testing_with_prefix.snakefile"
    prefix:
        'cv_split{cv_split}/deeprvat'
    config:
        config['deeprvat']

use rule * from deeprvat_workflow exclude config, choose_training_genes, train_bagging, regress, compute_burdens, compute_burdens_test, best_bagging_run, cleanup_burden_cache, link_burdens, link_burdens_test, all_burdens  as deeprvat_*

use rule evaluate from deeprvat_workflow as deeprvat_evaluate with:
    params:
        prefix = 'cv_split{cv_split}/deeprvat',
        use_seed_genes = '--use-seed-genes'

use rule regress from deeprvat_workflow as deeprvat_regress with:
    input:
        config = "cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/hpopt_config.yaml",
        chunks = lambda wildcards: expand(
            ('cv_split{{cv_split}}/deeprvat/{{phenotype}}/deeprvat/burdens/chunk{chunk}.' +
             ("finished" if wildcards.phenotype == burden_phenotype else "linked")),
            chunk=range(n_burden_chunks)
        ),
        phenotype_0_chunks =  expand(
            'cv_split{{cv_split}}/deeprvat/' + burden_phenotype + '/deeprvat/burdens/chunk{chunk}.finished',
            chunk=range(n_burden_chunks)
        ),
    params:
        prefix = 'cv_split{cv_split}/deeprvat'



use rule link_burdens from deeprvat_workflow as deeprvat_link_burdens with:
    params:
        prefix = 'cv_split{cv_split}/deeprvat'

use rule compute_burdens from deeprvat_workflow as deeprvat_compute_burdens with:
    params:
        prefix = 'cv_split{cv_split}/deeprvat'


rule all_training:
    input:
        expand('cv_split{cv_split}/deeprvat/models/repeat_{repeat}/best/bag_{bag}.ckpt',
               bag=range(n_bags), repeat=range(n_repeats),
               cv_split = range(cv_splits)),
        expand('cv_split{cv_split}/deeprvat/models/repeat_{repeat}/config.yaml',
               repeat=range(n_repeats),
               cv_split = range(cv_splits))

use rule best_training_run from deeprvat_workflow as deeprvat_best_training_run with:
    params:
        prefix = 'cv_split{cv_split}/deeprvat'

use rule train from deeprvat_workflow as deeprvat_train with:
    params:
        prefix = 'cv_split{cv_split}/deeprvat',
        phenotypes = " ".join( #TODO like need the prefix here as well
            [f"--phenotype {p} "
             f"cv_split{{cv_split}}/deeprvat/{p}/deeprvat/input_tensor.zarr "
             f"cv_split{{cv_split}}/deeprvat/{p}/deeprvat/covariates.zarr "
             f"cv_split{{cv_split}}/deeprvat/{p}/deeprvat/y.zarr"
             for p in phenotypes])

use rule choose_training_genes from deeprvat_workflow as deeprvat_choose_training_genes with:
    params:
        prefix = 'cv_split{cv_split}/deeprvat'

use rule config from deeprvat_workflow as deeprvat_config with:
    input:
        config = 'cv_split{cv_split}/deeprvat/config.yaml', # TODO: change this into cv specific config
        baseline = 'cv_split{cv_split}/baseline/{phenotype}/eval/burden_associations.parquet',
    params:
        baseline_results = '--baseline-results cv_split{cv_split}/baseline/{phenotype}/eval/burden_associations.parquet '


############################### Computation of test set deeprvat burdens ##############################################################
############################################################################################################################
rule make_deeprvat_test_config:
    input:
        config_train = 'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/hpopt_config.yaml'
    output:
        config_test = 'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/hpopt_config_test.yaml'
    shell:
        ' && '.join([
            conda_check,
            py_pipeline + 'cv_utils.py generate-test-config '
            '--fold {wildcards.cv_split} '
            '-m deeprvat '
            ' {input.config_train} {output.config_test}'
        ])

use rule association_dataset from deeprvat_workflow as deeprvat_association_dataset_test with:
    input:
        config = 'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/hpopt_config_test.yaml'
    output:
        'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/association_dataset_test.pkl'
    threads: 4



rule link_burdens_test:
    priority: 1
    input:
        checkpoints = lambda wildcards: [
            f'cv_split{{cv_split}}/deeprvat/models/repeat_{repeat}/best/bag_{bag}.ckpt'
            for repeat in range(n_repeats) for bag in range(n_bags)
        ],
        dataset = 'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/association_dataset_test.pkl',
        data_config = 'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/hpopt_config.yaml',
        model_config = 'cv_split{cv_split}/deeprvat/models/repeat_0/config.yaml',
    output:
        'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/burdens_test/chunk{chunk}.linked'
    threads: 8
    resources:
        mem_mb = lambda wildcards: 16384,
        load = lambda wildcards: 16000,
    shell:
        ' && '.join([
            ('deeprvat_associate compute-burdens '
             + debug +
             ' --n-chunks '+ str(n_burden_chunks) + ' '
             f'--link-burdens ../../../{burden_phenotype}/deeprvat/burdens_test/burdens.zarr '
             '--chunk {wildcards.chunk} '
             '--dataset-file {input.dataset} '
             '{input.data_config} '
             '{input.model_config} '
             '{input.checkpoints} '
             'cv_split{wildcards.cv_split}/deeprvat/{wildcards.phenotype}/deeprvat/burdens_test'),
            'touch {output}'
        ])

rule compute_burdens_test:
    priority: 10
    input:
        checkpoints = lambda wildcards: [
            f'cv_split{{cv_split}}/deeprvat/models/repeat_{repeat}/best/bag_{bag}.ckpt'
            for repeat in range(n_repeats) for bag in range(n_bags)
        ],
        dataset = 'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/association_dataset_test.pkl',
        data_config = 'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/hpopt_config.yaml',
        model_config = 'cv_split{cv_split}/deeprvat/models/repeat_0/config.yaml',
    output:
        'cv_split{cv_split}/deeprvat/{phenotype}/deeprvat/burdens_test/chunk{chunk}.finished'
    threads: 8
    resources:
        mem_mb = 2000000,        # Using this value will tell our modified lsf.profile not to set a memory resource
        load = 8000,
        gpus = 1
    shell:
        ' && '.join([
            ('deeprvat_associate compute-burdens '
             + debug +
             ' --n-chunks '+ str(n_burden_chunks) + ' '
             '--chunk {wildcards.chunk} '
             '--dataset-file {input.dataset} '
             '{input.data_config} '
             '{input.model_config} '
             '{input.checkpoints} '
             'cv_split{wildcards.cv_split}/deeprvat/{wildcards.phenotype}/deeprvat/burdens_test'),
            'touch {output}'
        ])

#dummy rule to indicate that all burdens.zarr files are there/complete for phenotype_prediction_pipeline.smk
rule all_burdens_zarr:
    input:
        expand('cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens_test/chunk{c}.finished',
                p = burden_phenotype, c = range(n_burden_chunks), cv_split=range(cv_splits)),
        expand('cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens_test/chunk{c}.linked',
                p = [pheno for pheno in phenotypes_testing if pheno != burden_phenotype],
                 c = range(n_burden_chunks), cv_split=range(cv_splits)),        
        expand('cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens/chunk{c}.finished',
                p = burden_phenotype, c = range(n_burden_chunks), cv_split=range(cv_splits)),
        expand('cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens/chunk{c}.linked',
                p = [pheno for pheno in phenotypes_testing if pheno != burden_phenotype],
                 c = range(n_burden_chunks), cv_split=range(cv_splits)),
    output:
        expand('cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens_test/burdens.zarr',
                p = burden_phenotype, cv_split=range(cv_splits)),
        expand('cv_split{cv_split}/deeprvat/{p}/deeprvat/burdens/burdens.zarr',
                p = burden_phenotype, cv_split=range(cv_splits)),
   

rule all_alternative_burdens_zarr:
    input:
        expand('cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/burdens/linking.finished',
            cv_split=range(cv_splits),
            btype = btypes,
            phenotype = phenotypes_testing),
    output:
        expand('cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/burdens/burdens.zarr',
            cv_split=range(cv_splits),
            btype = btypes,
            phenotype = phenotypes_testing),      
        expand('cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/burdens_test/burdens.zarr',
                cv_split=range(cv_splits),
                btype = btypes,
                phenotype = phenotypes_testing),     
#
############################### Computation of alternative burdens ##############################################################
############################################################################################################################


rule make_alt_burdens_test_config:
    input:
        config_train = 'cv_split{cv_split}/alternative_burdens/config.yaml'
    output:
        config_test = 'cv_split{cv_split}/alternative_burdens/config_test.yaml'
    shell:
        ' && '.join([
            conda_check,
            py_pipeline + 'cv_utils.py generate-test-config '
            '--fold {wildcards.cv_split} '
            '-m alternative_burdens '
            ' {input.config_train} {output.config_test}'
        ])

rule alt_burdens_association_dataset:
    input:
        config = lambda wildcards: 'cv_split{cv_split}/alternative_burdens/config.yaml' if wildcards.test_suffix == 'train' else 'cv_split{cv_split}/alternative_burdens/config_{test_suffix}.yaml'
    output:
        'cv_split{cv_split}/alternative_burdens/{btype}/association_dataset_{test_suffix}.pkl'
    threads: 4
    priority: 12
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * (attempt + 1),
        load = 64000
    shell:
        ' && '.join([
            conda_check,
            (py_pipeline + 'alternative_burdens.py make-dataset '
            '--data-key alt_burdens_data '
            '--btype {wildcards.btype} '
             + debug +
             '{input.config} '
             '{output}')
        ])

rule compute_alt_burdens_burdens:
    priority: 10
    input:
        config = lambda wildcards: 'cv_split{cv_split}/alternative_burdens/config.yaml' if wildcards.test_suffix == 'train' else 'cv_split{cv_split}/alternative_burdens/config_{test_suffix}.yaml',
        dataset = 'cv_split{cv_split}/alternative_burdens/{btype}/association_dataset_{test_suffix}.pkl'
    output:
        'cv_split{cv_split}/alternative_burdens/{btype}/burdens_{test_suffix}/chunk{altchunk}.finished'
    params:
        prefix= '.'
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 + (attempt - 1) * 2 * 8098,     
        load = 8000,
        # gpus = 1
    shell:
        ' && '.join([
            conda_check,
            cuda_visible_devices,
            (py_pipeline + 'alternative_burdens.py compute-alternative-burdens '
             + debug +
             ' --n-chunks '+ str(alt_burdens_chunks) + ' '
             '--chunk {wildcards.altchunk} '
             '--dataset-file {input.dataset} '
             '--btype {wildcards.btype} '
             '{input.config} '
             '{params.prefix}/cv_split{wildcards.cv_split}/alternative_burdens/{wildcards.btype}/burdens_{wildcards.test_suffix}'),
            'touch {output}'
        ])


rule move_alt_burdens_train_burdens:
    input:
        expand('cv_split{{cv_split}}/alternative_burdens/{{btype}}/burdens_train/chunk{altchunk}.finished',
                altchunk = range(alt_burdens_chunks))
    output:
       'cv_split{cv_split}/alternative_burdens/{btype}/burdens_train_moved.finished'
    shell:
        ' && '.join([
            'mv cv_split{wildcards.cv_split}/alternative_burdens/{wildcards.btype}/burdens_train/ cv_split{wildcards.cv_split}/alternative_burdens/{wildcards.btype}/burdens ',
            'touch {output}'
        ])

rule link_all_alt_burdens_burdens:
    input:
        expand('cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/burdens/linking.finished',
                cv_split=range(cv_splits),
                altchunk = range(alt_burdens_chunks),
                btype = btypes,
                phenotype = phenotypes_testing),
 #uses x and y from deeprvat directory because there phentoype specific datasets have been created
rule link_alt_burdens_burdens:
    input:
        test_burdens = expand('cv_split{{cv_split}}/alternative_burdens/{{btype}}/burdens_test/chunk{altchunk}.finished',
                altchunk = range(alt_burdens_chunks)),
        train_burdens = 'cv_split{cv_split}/alternative_burdens/{btype}/burdens_train_moved.finished'
    output:
       finished_train = 'cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/burdens/linking.finished',
       finished_test = 'cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/burdens_test/linking.finished',      
    params:
        link_path = 'cv_split{cv_split}/alternative_burdens/{phenotype}/{btype}/'
    shell:
        ' && '.join([
            'rm -rf {params.link_path} ',
            'mkdir -p {params.link_path}/burdens ',
            'mkdir -p {params.link_path}/burdens_test ',
            'ln -rsf cv_split{wildcards.cv_split}/deeprvat/{wildcards.phenotype}/deeprvat/burdens/x.zarr {params.link_path}/burdens/x.zarr ',
            'ln -rsf cv_split{wildcards.cv_split}/deeprvat/{wildcards.phenotype}/deeprvat/burdens/y.zarr {params.link_path}/burdens/y.zarr ',
            'ln -rsf cv_split{wildcards.cv_split}/alternative_burdens/{wildcards.btype}/burdens/burdens.zarr {params.link_path}/burdens/burdens.zarr ',
            'ln -rsf cv_split{wildcards.cv_split}/alternative_burdens/{wildcards.btype}/burdens/genes.npy {params.link_path}/burdens/genes.npy ',
            
            'ln -rsf cv_split{wildcards.cv_split}/deeprvat/{wildcards.phenotype}/deeprvat/burdens_test/x.zarr {params.link_path}/burdens_test/x.zarr ',
            'ln -rsf cv_split{wildcards.cv_split}/deeprvat/{wildcards.phenotype}/deeprvat/burdens_test/y.zarr {params.link_path}/burdens_test/y.zarr ',
            'ln -rsf cv_split{wildcards.cv_split}/alternative_burdens/{wildcards.btype}/burdens_test/burdens.zarr {params.link_path}/burdens_test/burdens.zarr ',
            'ln -rsf cv_split{wildcards.cv_split}/alternative_burdens/{wildcards.btype}/burdens_test/genes.npy {params.link_path}/burdens_test/genes.npy ',

            'touch {output.finished_train}' ,
            'touch {output.finished_test}' 
        ])


