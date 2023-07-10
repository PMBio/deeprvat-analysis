BTYPES_DICT = {'is_plof': 'plof',
              'CADD_raw':'cadd',
              'AbSplice_DNA': 'absplice',
              'PrimateAI_score': 'primateai',
              'sift_score': 'sift',
              'polyphen_score': 'polyphen',
              'SpliceAI_delta_score': 'splicai',
              'Consequence_missense_variant': 'missense'}


PLOF_CONSEQUENCES = [f'Consequence_{c}'  for c in ('splice_acceptor_variant', 'splice_donor_variant',
            'frameshift_variant', 'stop_gained', 'stop_lost', 'start_lost')]


PHENOTYPE_MAPPING = {'Forced_expiratory_volume_in_1_second_FEV1': 'FEV1',
                   'Mean_platelet_thrombocyte_volume': 'MPTVS',
                   'Body_mass_index_BMI': 'BMI',
                   'IGF_1': 'IGF-1',
                   'Red_blood_cell_erythrocyte_count': 'Erythrocyte count'}