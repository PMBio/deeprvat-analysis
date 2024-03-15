from pathlib import Path
from typing import List

py="python [path_to_deeprvat_analysis]/feature_importance/mutagenesis/"

path_to_inputs="[path_to_example_input_dir]"
pretrained_dir="[path_to_deeprvat]/pretrained_models"
config_file="[path_to_deeprvat_analysis]/feature_importance/config.yaml"


repeats = ['0', '1', '2', '3', '4', '5']
annotations = ['4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14']
# the numbers indicate 0-indexed order of the annotations used in training (config_file)  

 
rule all:
    input:
        expand('annot_{annot}/diff_repeat={repeat}.pkl', 
                annot=annotations, repeat=repeats),


rule run_over_bags:
    input:
        config_file = ''.join([config_file ]),
        input_dir = ''.join([ path_to_inputs ]),
        checkpoint_files = ''.join([ pretrained_dir + '/repeat_{repeat}/best']),  
    output:
        'annot_{annot}/diff_repeat={repeat}.pkl',
    threads: 1
    resources:
        mem_mb = 50000,
        disk_mb = 50000,
    shell:
        ' && '.join([
            py + 'explain_mutagenesis.py mutagenesis ' +
            '{input.config_file} ' +
            '{input.checkpoint_files} ' +
            '{input.input_dir} ' +
            'repeat={wildcards.repeat} ' +
            '{wildcards.annot} ' 
        ])

