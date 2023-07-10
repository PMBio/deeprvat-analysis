from pathlib import Path
from typing import List

configfile: "config.yaml"

# where explain_mutagenesis_indv.py file is saved along with pl_models.py
py = 'python ~/genopheno/genopheno/aggregation_metrics/'

# base_dir is where model ckpts are saved 
base_dir = '~/experiments/rvat/multipheno_bagging_reverse'

repeats = ['0', '1', '2', '3', '4', '5']
annotations = ['4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14']
   
rule all:
    input:
        expand('annot_{annot}/diff_repeat={repeat}.pkl', 
                annot=annotations, repeat=repeats),


rule run_over_bags:
    input:
        config_file = ''.join([ base_dir + '/models/repeat_{repeat}/config.yaml']),
        input_dir = ''.join([ base_dir ]),
        checkpoint_files = ''.join([ base_dir + '/models/repeat_{repeat}/best']),  
    output:
        'annot_{annot}/diff_repeat={repeat}.pkl',
    threads: 1
    resources:
        mem_mb = 50000,
        disk_mb = 50000,
    shell:
        ' && '.join([
            py + 'explain_mutagenesis_indv.py mutagenesis ' +
            '{input.config_file} ' +
            '{input.checkpoint_files} ' +
            '{input.input_dir} ' +
            'repeat={wildcards.repeat} ' +
            '{wildcards.annot} ' 
        ])

