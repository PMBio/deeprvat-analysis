import numpy as np
import pandas as pd
from plotnine import *
import glob
import yaml
import pickle
import os

SAMPLING_COLORS_LIST_10 = ['#c6dbef','#9ecae1', '#99CCCC', '#6baed6','#3182bd',
                     '#08519c',  '#8B4513', '#CD853F', '#D2691E', '#DEB887']

ANNO_COLOR_INDV_LIST =  ['#3182bd', '#6baed6',
                     '#b2723f', '#df880f',  '#8B4513', '#CD853F', '#D2691E', '#DEB887',
                     '#99CCCC', '#99CCCC', '#99CCCC', '#99CCCC', '#99CCCC', '#99CCCC', 
                     '#08519c',  '#08519c', '#08519c', '#08519c', '#08519c', '#08519c',
                   
                     '#c68eb7', '#c68eb7',
                     '#8370ab' , '#8370ab' 
                      ]


QUANT_COLOR_LIST =  ['#3182bd', '#6baed6',
                     '#b2723f', '#df880f',  '#8B4513', '#CD853F', '#D2691E', '#DEB887',
                     '#99CCCC', '#99CCCC', '#99CCCC', '#99CCCC', '#08519c', '#08519c', 
                     '#08519c',  '#08519c', '#08519c', '#08519c', '#08519c', '#08519c',
                   
                     '#c68eb7', '#c68eb7',
                     #'#8370ab' , '#8370ab' 
                      ]


#'#3182bd',  MAF
#'#6baed6', CADD
# '#99CCCC', protein function 
# 08519c pLof
#'#b2723f', BHLHE
#'#df880f', CTCF
#'#8B4513', ccFos
#'#CD853F', H3k4me2, H3kme9
#'#D2691E', H3k4me2
#'#DEB887', p300
# #c68eb7  inframe indels
# '#8370ab'  splicing 


ANNOTATION_NAMES = ['combined_UKB_NFE_AF_MB', 'CADD_raw',
         'sift_score', 'polyphen_score',
         'Consequence_splice_acceptor_variant', 'Consequence_splice_donor_variant',
         'Consequence_stop_gained', 'Consequence_frameshift_variant',
         'Consequence_stop_lost','Consequence_start_lost',
         'Consequence_inframe_insertion','Consequence_inframe_deletion',
         'Consequence_missense_variant', 'Consequence_protein_altering_variant',
         'Consequence_splice_region_variant', 'condel_score',
         'DeepSEA_PC_1','DeepSEA_PC_2', 'DeepSEA_PC_3', 'DeepSEA_PC_4',
         'DeepSEA_PC_5', 'DeepSEA_PC_6','PrimateAI_score', 'AbSplice_DNA',
         'DeepRipe_plus_QKI_lip_hg2', 'DeepRipe_plus_QKI_clip_k5',
         'DeepRipe_plus_KHDRBS1_clip_k5', 'DeepRipe_plus_ELAVL1_parclip',
         'DeepRipe_plus_TARDBP_parclip','DeepRipe_plus_HNRNPD_parclip',
         'DeepRipe_plus_MBNL1_parclip', 'DeepRipe_plus_QKI_parclip',
         'SpliceAI_delta_score']

BINARY_ANNOTATION_CATEGORIES = { 
                          'Consequence_missense_variant': 'Protein function',
                          'Consequence_protein_altering_variant': 'Protein function',
                          'Consequence_stop_lost': 'pLof', 
                          'Consequence_start_lost': 'pLof', 
                          'Consequence_splice_acceptor_variant': 'pLof',
                          'Consequence_splice_donor_variant': 'pLof', 
                          'Consequence_stop_gained': 'pLof',
                          'Consequence_frameshift_variant': 'pLof',
                          'Consequence_splice_region_variant': 'Splicing',
                          'Consequence_inframe_insertion': 'Inframe indels',
                          'Consequence_inframe_deletion': 'Inframe indels',

                        }


ANNOTATION_CATEGORIES = { 'combined_UKB_NFE_AF_MB': 'UKB MAF',
                          'CADD_raw': 'CADD raw',
                          'DeepSEA_PC_1': 'BHLHE40',
                          'DeepSEA_PC_2': 'CTCF',
                          'DeepSEA_PC_3': 'cFos',
                          'DeepSEA_PC_4': 'H3K4me2, H3k9ac',
                          'DeepSEA_PC_5': 'H3K4me2',
                          'DeepSEA_PC_6': 'p300',
                         
                          #'Consequence_missense_variant': 'Protein function',
                          'sift_score': 'Protein function',
                          'polyphen_score' : 'Protein function', 
                          'PrimateAI_score': 'Protein function',
                          'condel_score' : 'Protein function', 
                          #'Consequence_protein_altering_variant': 'Protein function',
                         
                          'DeepRipe_plus_QKI_lip_hg2': 'RNA-binding',
                          'DeepRipe_plus_QKI_clip_k5': 'RNA-binding',
                          'DeepRipe_plus_KHDRBS1_clip_k5': 'RNA-binding',
                          'DeepRipe_plus_ELAVL1_parclip': 'RNA-binding',
                          'DeepRipe_plus_TARDBP_parclip': 'RNA-binding',
                          'DeepRipe_plus_HNRNPD_parclip': 'RNA-binding',
                          'DeepRipe_plus_MBNL1_parclip': 'RNA-binding',
                          'DeepRipe_plus_QKI_parclip': 'RNA-binding',
                         
                          #'Consequence_stop_lost': 'pLof', 
                          #'Consequence_start_lost': 'pLof', 
                          #'Consequence_splice_acceptor_variant': 'pLof',
                          #'Consequence_splice_donor_variant': 'pLof', 
                          #'Consequence_stop_gained': 'pLof',
                          #'Consequence_frameshift_variant': 'pLof',
                         
                          #'Consequence_splice_region_variant': 'Splicing',
                          'AbSplice_DNA': 'Splicing',
                          'SpliceAI_delta_score': 'Splicing',
                         
                          #'Consequence_inframe_insertion': 'Inframe indels',
                          #'Consequence_inframe_deletion': 'Inframe indels',

                        }





ANNOTATION_CODES = {      'combined_UKB_NFE_AF_MB': 1,
                          'CADD_raw': 2,
                          'DeepSEA_PC_1': 3,
                          'DeepSEA_PC_2': 4,
                          'DeepSEA_PC_3': 5,
                          'DeepSEA_PC_4': 6,
                          'DeepSEA_PC_5': 7,
                          'DeepSEA_PC_6': 8,
                    
                          #'Consequence_missense_variant': 9,
                          'sift_score': 9,
                          'polyphen_score' : 9, 
                          'PrimateAI_score': 9,
                          'condel_score' : 9, 
                          #'Consequence_protein_altering_variant': 9,
                    
                          'DeepRipe_plus_QKI_lip_hg2': 10,
                          'DeepRipe_plus_QKI_clip_k5': 10,
                          'DeepRipe_plus_KHDRBS1_clip_k5': 10,
                          'DeepRipe_plus_ELAVL1_parclip':10,
                          'DeepRipe_plus_TARDBP_parclip': 10,
                          'DeepRipe_plus_HNRNPD_parclip': 10,
                          'DeepRipe_plus_MBNL1_parclip': 10,
                          'DeepRipe_plus_QKI_parclip': 10,
                    
                          #'Consequence_stop_lost': 11, 
                          #'Consequence_start_lost': 11, 
                          #'Consequence_splice_acceptor_variant': 11,
                          #'Consequence_splice_donor_variant': 11, 
                          #'Consequence_stop_gained': 11,
                          #'Consequence_frameshift_variant': 11,
                    
                          'AbSplice_DNA': 12,
                          #'Consequence_splice_region_variant': 12,
                          'SpliceAI_delta_score': 12, 
                          #'Consequence_inframe_insertion': 13,
                          #'Consequence_inframe_deletion': 13

                        }




ANNOTATION_PRINT = { 'combined_UKB_NFE_AF_MB': 'UKB MAF',
                          'CADD_raw': 'CADD raw',
                          'DeepSEA_PC_1': 'BHLHE40',
                          'DeepSEA_PC_2': 'CTCF',
                          'DeepSEA_PC_3': 'cFos',
                          'DeepSEA_PC_4': 'H3K4me2, H3k9ac',
                          'DeepSEA_PC_5': 'H3K4me2',
                          'DeepSEA_PC_6': 'p300',
                    
                          #'Consequence_missense_variant': 'missense variant',
                    
                          'sift_score': 'sift score',
                          'polyphen_score' : 'polyphen score', 
                          'PrimateAI_score': 'PrimateAI score',
                          'condel_score' : 'condel score', 
                    
                          #'Consequence_protein_altering_variant': 'protein altering variant',
                          #'Consequence_stop_lost': 'stop lost', 
                          #'Consequence_start_lost': 'start lost', 
                          #'Consequence_splice_acceptor_variant': 'splice acceptor variant',
                          #'Consequence_splice_donor_variant': 'splice donor variant', 
                          #'Consequence_stop_gained': 'stop gained',
                          #'Consequence_frameshift_variant': 'frameshift variant',
                    
                          'DeepRipe_plus_QKI_lip_hg2': 'RBP-QKI1',
                          'DeepRipe_plus_QKI_clip_k5': 'RBP-QKI2',
                          'DeepRipe_plus_KHDRBS1_clip_k5': 'RBP-KHDRBS1', 
                          'DeepRipe_plus_ELAVL1_parclip': 'RBP-ELAVL1',
                          'DeepRipe_plus_TARDBP_parclip': 'RBP-TARDBP', 
                          'DeepRipe_plus_HNRNPD_parclip': 'RBP-HNRNPD', 
                          'DeepRipe_plus_MBNL1_parclip': 'RBP-MBNL1', 
                          'DeepRipe_plus_QKI_parclip': 'RBP-QKI3', 
                        
                    
                          #'Consequence_splice_region_variant': 'splice region variant',
                          'AbSplice_DNA': 'AbSplice DNA',
                          'SpliceAI_delta_score': 'SpliceAI',
                    
                          #'Consequence_inframe_insertion': 'inframe insertion',
                          #'Consequence_inframe_deletion': 'inframe deletion',

                        }



    

def update_annots_importance_order(df, 
                                   col_title='anno_category',
                                   val_title='value'):
    agg_annots = df.groupby(col_title)[val_title].agg(['mean'])
    sorted_annots = list(agg_annots.sort_values(by='mean', ascending=False).index)
    df[col_title] = df[col_title].astype('category')
    df[col_title] = df[col_title].cat.reorder_categories(sorted_annots)
    
    return df





def collect_repeats_in_one_sampling(sampling_dir, 
                       annotation_names, 
                       repeats=6):
    importance_by_repeat = {f'repeat_{i}':{} for i in range(repeats)}
    sampling_no = os.path.basename(sampling_dir[:-1])
    print(f'Collecting per repeat/per pheno importance: {sampling_no}')
    for repeat_pth in glob.glob(f'{sampling_dir}*annots.pkl'):
        repeat_no = os.path.basename(repeat_pth)[:8]
        with open(repeat_pth, "rb") as f:
            repeat_dict = pickle.load(f)

        per_pheno_importance = {p:[] for p in repeat_dict.keys()}
        for pheno in repeat_dict.keys():
            ## using absolute values
            pheno_importance = np.absolute(repeat_dict[pheno])
            agg_over_variants = np.mean(pheno_importance, axis=3)
            agg_over_genes = np.mean(agg_over_variants, axis=1)
            agg_over_samples = np.mean(agg_over_genes, axis=0)
            # save per pheno annotation importance
            per_pheno_importance[pheno] = agg_over_samples
        # save per repeat per pheno importance
        importance_by_repeat[repeat_no] = per_pheno_importance
    
    return importance_by_repeat



def agg_over_repeats(per_repeat_importance):
    pheno_over_repeats = {p:[] 
                for p in per_repeat_importance['repeat_0'].keys()}
    pheno_agg_repeats = {}
    for repeat in per_repeat_importance.keys():
        per_pheno = per_repeat_importance[repeat]
        for pheno in per_pheno.keys():
            pheno_over_repeats[pheno].append(per_pheno[pheno])

    for pheno in pheno_over_repeats.keys():
        pheno_agg_repeats[pheno] = np.mean(pheno_over_repeats[pheno], axis=0)
        
    return pheno_agg_repeats


    
def add_plot_helper_columns(df):
    assert ('annotation' in df.columns) == True, 'missing annotation column'
    df['anno_code'] = df.apply(lambda x: 
                                       ANNOTATION_CODES[x['annotation']], axis=1)
    df['anno_category'] = df.apply(lambda x: 
                                       ANNOTATION_CATEGORIES[x['annotation']], axis=1)
    df['anno_print'] = df.apply(lambda x: 
                                        ANNOTATION_PRINT[x['annotation']], axis=1)
    
    return df 