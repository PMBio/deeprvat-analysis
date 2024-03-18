from pathlib import Path
from typing import List


py="python [path_to_deeprvat_analysis]/feature_importance/shap/"

path_to_inputs="[path_to_example_input_dir]"
pretrained_dir="[path_to_deeprvat]/pretrained_models"
config_file="[path_to_deeprvat_analysis]/feature_importance/config.yaml"



repeats = ['0', '1','2','3','4','5']
samplings = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15']

   
rule all:
    input:
        expand('sample_{sampling}/repeat_{repeat}_shap_avg_annots.pkl', 
                sampling=samplings, repeat=repeats),


    
rule run_over_bags:
    input:
        config_file = ''.join([ config_file ]),
        input_dir = ''.join([ path_to_inputs ]),
        checkpoint_files = ''.join([ pretrained_dir + '/repeat_{repeat}/best']),                                           
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

