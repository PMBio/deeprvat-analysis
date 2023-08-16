from pathlib import Path
from typing import List

configfile: 'config.yaml'

debug_flag = config.get('debug', False)
phenotypes = config['phenotypes']
phenotypes = phenotypes.keys() if type(phenotypes) == dict else phenotypes
burden_phenotype = phenotypes[0]
if type(phenotypes) == dict:
    training_phenotypes = {p: config['phenotypes'][p].get('training_phenotype', p)
                        for p in phenotypes}
else:
    training_phenotypes = phenotypes
n_burden_chunks = config.get('n_burden_chunks', 1) if not debug_flag else 2
n_regression_chunks = config.get('n_regression_chunks', 40) if not debug_flag else 2
n_trials = config['hyperparameter_optimization']['n_trials']
n_bags = config['training']['n_bags'] if not debug_flag else 3
n_repeats = config['n_repeats']
debug = '--debug ' if debug_flag else ''
do_scoretest = '--do-scoretest ' if config.get('do_scoretest', False) else ''
tensor_compression_level = config['training'].get('tensor_compression_level', 1)

#######################

#this pipeline is a modfication of training_association_testing.snakefile which adds prefixes to rules using 
#relative paths in their script to ensure that the pipeline can be used as a module from other snakefiles.
#######################


wildcard_constraints:
    repeat="\d+",
    trial="\d+",

rule all:
    input:
        expand("{phenotype}/deeprvat/eval/significant.parquet",
               phenotype=phenotypes),
        expand("{phenotype}/deeprvat/eval/all_results.parquet",
               phenotype=phenotypes)

rule evaluate:
    input:
        associations = expand('{{phenotype}}/deeprvat/repeat_{repeat}/results/burden_associations.parquet',
                              repeat=range(n_repeats)),
        config = '{phenotype}/deeprvat/hpopt_config.yaml',
    output:
        "{phenotype}/deeprvat/eval/significant.parquet",
        "{phenotype}/deeprvat/eval/all_results.parquet"
    threads: 1
    params:
        prefix = '.',
        use_seed_genes = '--use-seed-genes'
    resources:
        mem_mb = 4098
    shell:
        'deeprvat_evaluate '
        + debug +
        '{params.use_seed_genes} '
        # '--use-seed-genes '
        '--n-repeats {n_repeats} '
        '--correction-method FDR '
        '{input.associations} '
        '{input.config} '
        '{params.prefix}/{wildcards.phenotype}/deeprvat/eval'

rule all_regression:
    input:
        expand('{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations.parquet',
               phenotype=phenotypes, type=['deeprvat'], repeat=range(n_repeats)),

rule combine_regression_chunks:
    input:
        expand('{{phenotype}}/deeprvat/repeat_{{repeat}}/results/burden_associations_{chunk}.parquet', chunk=range(n_regression_chunks)),
    output:
        '{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations.parquet',
    threads: 1
    shell:
        'deeprvat_associate combine-regression-results '
        '--model-name repeat_{wildcards.repeat} '
        '{input} '
        '{output}'

rule regress:
    input:
        config = "{phenotype}/deeprvat/hpopt_config.yaml",
        chunks = lambda wildcards: expand(
            ('{{phenotype}}/deeprvat/burdens/chunk{chunk}.' +
             ("finished" if wildcards.phenotype == burden_phenotype else "linked")),
            chunk=range(n_burden_chunks)
        ),
        phenotype_0_chunks =  expand(
            burden_phenotype + '/deeprvat/burdens/chunk{chunk}.finished',
            chunk=range(n_burden_chunks)
        ),
    params:
        prefix = '.'
    output:
        temp('{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations_{chunk}.parquet'),
    resources:
        mem_mb = lambda wildcards, attempt: 24572 + (attempt - 1) * 4098,
        load = lambda wildcards, attempt: 24000 + (attempt - 1) * 4000
    threads: 2
    shell:
        'deeprvat_associate regress '
        + debug +
        '--chunk {wildcards.chunk} '
        '--n-chunks ' + str(n_regression_chunks) + ' '
        '--use-bias '
        '--repeat {wildcards.repeat} '
        + do_scoretest +
        '{input.config} '
        '{params.prefix}/{wildcards.phenotype}/deeprvat/burdens ' #TODO make this w/o repeats
        '{params.prefix}/{wildcards.phenotype}/deeprvat/repeat_{wildcards.repeat}/results'

rule all_burdens:
    input:
        [
            (f'{p}/deeprvat/burdens/chunk{c}.' +
             ("finished" if p == burden_phenotype else "linked"))
            for p in phenotypes
            for c in range(n_burden_chunks)
        ]

rule link_burdens:
    priority: 1
    input:
        checkpoints = lambda wildcards: [
            f'models/repeat_{repeat}/best/bag_{bag}.ckpt'
            for repeat in range(n_repeats) for bag in range(n_bags)
        ],
        dataset = '{phenotype}/deeprvat/association_dataset.pkl',
        data_config = '{phenotype}/deeprvat/hpopt_config.yaml',
        model_config = 'models/repeat_0/config.yaml',
    output:
        '{phenotype}/deeprvat/burdens/chunk{chunk}.linked'
    params:
        prefix = '.'
    resources:
        mem_mb = lambda wildcards: 16384,
        load = lambda wildcards: 16000,
    threads: 8
    shell:
        ' && '.join([
            ('deeprvat_associate compute-burdens '
             + debug +
             ' --n-chunks '+ str(n_burden_chunks) + ' '
             f'--link-burdens ../../../{burden_phenotype}/deeprvat/burdens/burdens.zarr '
             '--chunk {wildcards.chunk} '
             '--dataset-file {input.dataset} '
             '{input.data_config} '
             '{input.model_config} '
             '{input.checkpoints} '
             '{params.prefix}/{wildcards.phenotype}/deeprvat/burdens'),
            'touch {output}'
        ])

rule compute_burdens:
    priority: 10
    input:
        reversed = "models/reverse_finished.tmp",
        checkpoints = lambda wildcards: [
            f'models/repeat_{repeat}/best/bag_{bag}.ckpt'
            for repeat in range(n_repeats) for bag in range(n_bags)
        ],
        dataset = '{phenotype}/deeprvat/association_dataset.pkl',
        data_config = '{phenotype}/deeprvat/hpopt_config.yaml',
        model_config = 'models/repeat_0/config.yaml',
    output:
        '{phenotype}/deeprvat/burdens/chunk{chunk}.finished'
    params:
        prefix = '.'
    resources:
        mem_mb = 2000000,        # Using this value will tell our modified lsf.profile not to set a memory resource
        load = 8000,
        gpus = 1
    threads: 8
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
             '{params.prefix}/{wildcards.phenotype}/deeprvat/burdens'),
            'touch {output}'
        ])

rule all_association_dataset:
    input:
        expand('{phenotype}/deeprvat/association_dataset.pkl',
               phenotype=phenotypes)

rule association_dataset:
    input:
        config = '{phenotype}/deeprvat/hpopt_config.yaml'
    output:
        '{phenotype}/deeprvat/association_dataset.pkl'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * (attempt + 1),
        load = 64000
    shell:
        'deeprvat_associate make-dataset '
        + debug +
        '{input.config} '
        '{output}'

rule reverse_models:
    input:
        checkpoints = expand('models/repeat_{repeat}/best/bag_{bag}.ckpt',
                             bag=range(n_bags), repeat=range(n_repeats)),
        config = 'models/repeat_0/config.yaml/', 
    output:
        "models/reverse_finished.tmp"
    threads: 4
    shell:
        " && ".join([
            ("deeprvat_associate reverse-models "
             "{input.config} "
             "{input.checkpoints}"),
            "touch {output}"
        ])

rule all_training:
    input:
        expand('models/repeat_{repeat}/best/bag_{bag}.ckpt',
               bag=range(n_bags), repeat=range(n_repeats)),
        expand('models/repeat_{repeat}/config.yaml',
               repeat=range(n_repeats))

rule best_training_run:
    input:
        expand('models/repeat_{{repeat}}/trial{trial_number}/config.yaml',
               trial_number=range(n_trials)),
    output:
        checkpoints = expand('models/repeat_{{repeat}}/best/bag_{bag}.ckpt',
                             bag=range(n_bags)),
        config = 'models/repeat_{repeat}/config.yaml'
    params:
        prefix = '.'
    threads: 1
    shell:
        (
            'deeprvat_train best-training-run '
            + debug +
            '{params.prefix}/models/repeat_{wildcards.repeat} '
            '{params.prefix}/models/repeat_{wildcards.repeat}/best '
            '{params.prefix}/models/repeat_{wildcards.repeat}/hyperparameter_optimization.db '
            '{output.config}'
        )

rule train:
    input:
        config = expand('{phenotype}/deeprvat/hpopt_config.yaml',
                        phenotype=phenotypes),
        input_tensor = expand('{phenotype}/deeprvat/input_tensor.zarr',
                              phenotype=phenotypes),
        covariates = expand('{phenotype}/deeprvat/covariates.zarr',
                            phenotype=phenotypes),
        y = expand('{phenotype}/deeprvat/y.zarr',
                   phenotype=phenotypes),
    output:
        config = 'models/repeat_{repeat}/trial{trial_number}/config.yaml',
        finished = 'models/repeat_{repeat}/trial{trial_number}/finished.tmp'
    resources:
        mem_mb = 2000000,        # Using this value will tell our modified lsf.profile not to set a memory resource
        load = 8000,
        gpus = 1
    params:
        phenotypes = " ".join(
            [f"--phenotype {p} "
             f"{p}/deeprvat/input_tensor.zarr "
             f"{p}/deeprvat/covariates.zarr "
             f"{p}/deeprvat/y.zarr"
             for p in phenotypes]),
        prefix = '.'
    shell:
        ' && '.join([
            'deeprvat_train train '
            + debug +
            '--trial-id {wildcards.trial_number} '
            "{params.phenotypes} "
            'config.yaml '
            '{params.prefix}/models/repeat_{wildcards.repeat}/trial{wildcards.trial_number} '
            '{params.prefix}/models/repeat_{wildcards.repeat}/hyperparameter_optimization.db',
            'touch {output.finished}'
        ])

rule all_training_dataset:
    input:
        input_tensor = expand('{phenotype}/deeprvat/input_tensor.zarr',
                              phenotype=phenotypes, repeat=range(n_repeats)),
        covariates = expand('{phenotype}/deeprvat/covariates.zarr',
                            phenotype=phenotypes, repeat=range(n_repeats)),
        y = expand('{phenotype}/deeprvat/y.zarr',
                   phenotype=phenotypes, repeat=range(n_repeats))

rule training_dataset:
    input:
        config = '{phenotype}/deeprvat/hpopt_config.yaml',
        training_dataset = '{phenotype}/deeprvat/training_dataset.pkl'
    output:
        input_tensor = directory('{phenotype}/deeprvat/input_tensor.zarr'),
        covariates = directory('{phenotype}/deeprvat/covariates.zarr'),
        y = directory('{phenotype}/deeprvat/y.zarr')
    threads: 8
    priority: 50    
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * (attempt + 1),
        load = 16000
    shell:
        (
            'deeprvat_train make-dataset '
            + debug +
            '--compression-level ' + str(tensor_compression_level) + ' '
            '--training-dataset-file {input.training_dataset} '
            '{input.config} '
            '{output.input_tensor} '
            '{output.covariates} '
            '{output.y}'
        )

rule training_dataset_pickle:
    input:
        '{phenotype}/deeprvat/hpopt_config.yaml'
    output:
        '{phenotype}/deeprvat/training_dataset.pkl'
    threads: 1
    resources:
        mem_mb = 40000, # lambda wildcards, attempt: 38000 + 12000 * attempt
        load = 16000
    shell:
        (
            'deeprvat_train make-dataset '
            '--pickle-only '
            '--training-dataset-file {output} '
            '{input} '
            'dummy dummy dummy'
        )

rule all_config:
    input:
        seed_genes = expand('{phenotype}/deeprvat/seed_genes.parquet',
                            phenotype=phenotypes),
        config = expand('{phenotype}/deeprvat/hpopt_config.yaml',
                        phenotype=phenotypes),
        baseline = expand('{phenotype}/deeprvat/baseline_results.parquet',
                          phenotype=phenotypes),

rule config:
    input:
        config = 'config.yaml',
        baseline = lambda wildcards: [
            str(Path(r['base']) / wildcards.phenotype / r['type'] /
                'eval/burden_associations_testing.parquet')
            for r in config['baseline_results']
        ]
    output:
        seed_genes = '{phenotype}/deeprvat/seed_genes.parquet',
        config = '{phenotype}/deeprvat/hpopt_config.yaml',
        baseline = '{phenotype}/deeprvat/baseline_results.parquet',
    threads: 1
    params:
        baseline_results = lambda wildcards, input: ''.join([
            f'--baseline-results {b} '
            for b in input.baseline
        ])
    shell:
        (
            'deeprvat_config update-config '
            '--phenotype {wildcards.phenotype} '
            '{params.baseline_results}'
            '--baseline-results-out {output.baseline} '
            '--seed-genes-out {output.seed_genes} '
            '{input.config} '
            '{output.config}'
        )
