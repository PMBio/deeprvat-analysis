from pathlib import Path
from typing import List

configfile: 'config.yaml'

py = 'python ~/genopheno/genopheno/aggregation_metrics/'

## base_dir is where ckpts are saved
base_dir = '~/experiments/rvat/multipheno_bagging_reverse'
repeats = ['0', '1','2','3','4','5']
samplings = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15']

   
rule all:
    input:
        expand('sample_{sampling}/repeat_{repeat}_shap_avg_annots.pkl', 
                sampling=samplings, repeat=repeats),


    
rule run_over_bags:
    input:
        config_file = ''.join([ base_dir + '/models/repeat_{repeat}/config.yaml']),
        input_dir = ''.join([ base_dir ]),
        checkpoint_files = ''.join([ base_dir + '/models/repeat_{repeat}/best']),                                           
    output:
        'sample_{sampling}/repeat_{repeat}_shap_avg_annots.pkl'
    threads: 1
    resources:
        mem_mb = 100000,
        disk_mb = 100000,
        load = 600000
    shell:
        ' && '.join([
            py + 'explain_shap.py feature-importance ' +
            '{input.config_file} ' +
            '{input.checkpoint_files} ' +
            '{input.input_dir} ' +
            'repeat_{wildcards.repeat} ' +
            '{wildcards.sampling} ' 

        ])

