import numpy as np
import pandas as pd
from plotnine import *
import glob
import yaml
import pickle
import os


font_size = 8
font_family = "Helvetica"
    

ANNO_COLOR_INDV_LIST =  ['#3182bd', '#6baed6',
                     '#b2723f', '#df880f',  '#8B4513', '#CD853F', '#D2691E', '#DEB887',
                     '#99CCCC', '#99CCCC', '#99CCCC', '#99CCCC', '#99CCCC', '#99CCCC', 
                     '#08519c',  '#08519c', '#08519c', '#08519c', '#08519c', '#08519c',
                   
                     '#c68eb7', '#c68eb7',
                     '#8370ab' , '#8370ab' 
                      ]


BINARY_ANNO_COLOR_INDV_LIST = [

                     '#99CCCC', '#99CCCC', 
                     '#08519c',  '#08519c', '#08519c', '#08519c', '#08519c', '#08519c',
                     
                     '#8370ab', #splicing
                     '#c68eb7', '#c68eb7',
                     
                      ]

BINARY_ANNO_COLOR_GROUP_LIST = [
                                 '#99CCCC', 
                                 '#08519c',                      
                                 '#8370ab', #splicing
                                 '#c68eb7'                   
                                  ]

BINARY_ANNO_GROUP = ['Protein function', 'pLof', 'Splicing', 'Inframe indels']

BINARY_ANNOTATION_CODES = {      
                          'missense_variant': 9,
                          'protein_altering': 9,
                          'stop_lost': 10, 
                          'start_lost': 10, 
                          'splice_acceptor': 10,
                          'splice_donor': 10, 
                          'stop_gained': 10,
                          'frameshift': 10,
                          'splice_region': 11,
                          'inframe_insertion': 12,
                          'inframe_deletion': 12

                        }


ANNOTATION_CLEAN = {     'missense_variant': 'missense variant',
                          'protein_altering': 'protein altering variant',
                          'stop_lost': 'stop lost', 
                          'start_lost': 'start lost', 
                          'splice_acceptor': 'splice acceptor variant',
                          'splice_donor': 'splice donor variant', 
                          'stop_gained': 'stop gained',
                          'frameshift': 'frameshift variant',
                          'splice_region': 'splice region variant',
                          'inframe_insertion': 'inframe insertion',
                          'inframe_deletion': 'inframe deletion',

                        }



BINARY_ANNOTATIONS = {'annot_4': 'splice_acceptor', 
               'annot_5': 'splice_donor', 
               'annot_6': 'stop_gained', 
               'annot_7': 'frameshift', 
               'annot_8': 'stop_lost', 
               'annot_9': 'start_lost', 
               'annot_10': 'inframe_insertion', 
               'annot_11': 'inframe_deletion', 
               'annot_12': 'missense_variant', 
               'annot_13': 'protein_altering',
               'annot_14': 'splice_region'}




BINARY_ANNOTATION_CATEGORIES = { 
                          'missense_variant': 'Protein function',
                          'protein_altering': 'Protein function',
                          'stop_lost': 'pLof', 
                          'start_lost': 'pLof', 
                          'splice_acceptor': 'pLof',
                          'splice_donor': 'pLof', 
                          'stop_gained': 'pLof',
                          'frameshift': 'pLof',
                          'splice_region': 'Splicing',
                          'inframe_insertion': 'Inframe indels',
                          'inframe_deletion': 'Inframe indels',

                        }

PHENO_CODES = {'Apolipoprotein_A':1,
             'Apolipoprotein_B':1,
             'Calcium':1,
             'Cholesterol':1,
             'HDL_cholesterol':1,
             'IGF_1':1,
             'LDL_direct':1,
             'Lymphocyte_percentage':1,
             'MPTV':1,
             'Mean_corpuscular_volume':1,
             'Mean_reticulocyte_volume':1,
             'Neutrophill_count':2,
             'Platelet_count':2,
             'Platelet_crit':2,
             'Platelet_distribution_width':2,
             'Red_blood_cell_erythrocyte_count':2,
             'SHBG':2,
             'Standing_height':2,
             'Total_bilirubin':2,
             'Triglycerides':2,
             'Urate':2}

PHENO_NAMES = {'Apolipoprotein_A': 'Apolipo\nprotein\nA',
             'Apolipoprotein_B': 'Apolipo\nprotein\nB',
             'Calcium': 'Calcium',
             'Cholesterol': 'Cholesterol',
             'HDL_cholesterol': 'HDL',
             'IGF_1': 'IGF-1',
             'LDL_direct': 'LDL',
             'Lymphocyte_percentage': 'Lymphocyte\npercentage',
             'MPTV': 'MPTVS',
             'Mean_corpuscular_volume': 'Mean\ncorpuscular\nvolume',
             'Mean_reticulocyte_volume': 'Mean\nreticulocyte\nvolume',
             'Neutrophill_count': 'Neutrophill\ncount',
             'Platelet_count': 'Platelet\ncount',
             'Platelet_crit':'Platelet\ncrit',
             'Platelet_distribution_width':'Platelet\ndistr.\nwidth',
             'Red_blood_cell_erythrocyte_count':'Erythrocyte\ncount',
             'SHBG':'SHBG',
             'Standing_height':'Standing\nheight',
             'Total_bilirubin': 'Total\nbilirubin',
             'Triglycerides': 'Triglycerides',
             'Urate': 'Urate'}

plot_font = element_text(size = 8, family = "Helvetica")


### to observe variance across phenotypes
def save_plot_variance_across_pheno(df_relative_melted ):    
    p = (
            ggplot(df_relative_melted)
            + geom_boxplot(aes(x='reorder(anno_print, anno_code)', y='value', 
                           fill='variable')) ## variable # fill='variable'
            + theme_classic()
            + theme(axis_text_x=element_text(rotation=45, hjust=1),
                    #text=element_text(size=20),
                    figure_size=(10, 8),
                    text = plot_font,
                    axis_text = plot_font,
                    axis_title = plot_font,
                    legend_title=element_blank(),
                    legend_position='none'
                   )
            + scale_fill_manual(values=BINARY_ANNO_COLOR_INDV_LIST, 
                               breaks = list(BINARY_ANNOTATION_CODES.keys()),
                               labels = list(BINARY_ANNOTATION_CATEGORIES.values()),
                               #guide=True
                               )
            + labs(x='Annotation', y='Relative |absolute difference|')
        )
    ggsave(plot=p, filename='binary_importance_variance_across_pheno.pdf', 
                   limitsize=False, verbose = False)
    
    
    
def plot_individual_annot(df_relative_melted):
    plot_font = element_text(size = 8, family = "Helvetica")
    p = (
            ggplot(df_relative_melted, 
                   aes(x='reorder(anno_print, anno_code)', y='relative_importance', 
                           fill='annotation'))
            + geom_bar(stat = 'identity')
            + theme_classic()
            + theme(axis_text_x=element_text(rotation=45, hjust=1),
                    #text=element_text(size=20),
                    figure_size=(10, 8),
                    text = plot_font,
                    axis_text = plot_font,
                    axis_title = plot_font,
                    legend_title=element_blank(),
                    legend_position='top'
                   )
            + scale_fill_manual(values=BINARY_ANNO_COLOR_INDV_LIST, 
                               breaks = list(BINARY_ANNOTATION_CODES.keys()),
                               labels = list(BINARY_ANNOTATION_CATEGORIES.values()),
                               #guide=True
                               )
            + labs(x='Binary annotation', y='Relative |absolute difference|')
        )
    ggsave(plot=p, filename='binary_importance_by_annotation.pdf', 
                   limitsize=False, verbose = False)
    
    
    

    
    

    
def plot_category_annot(df_relative_melted):
    p = (
            ggplot(df_relative_melted, 
                   aes(x='reorder(category, anno_code)', y='relative_importance', 
                           fill='anno_code'))
            + geom_bar(stat = 'identity')
            + theme_classic()
            + theme(axis_text_x=element_text(rotation=45, hjust=1),
                    #text=element_text(size=20),
                    figure_size=(10, 8),
                    text = plot_font,
                    axis_text = plot_font,
                    axis_title = plot_font,
                    legend_title=element_blank(),
                    legend_position='none'
                   )
            + scale_fill_manual(values=BINARY_ANNO_COLOR_GROUP_LIST, 
                               breaks = BINARY_ANNO_GROUP,
                               labels = BINARY_ANNO_GROUP,
                               #guide=True
                               )
            + labs(x='Binary annotation groups',
                   y='Relative |average absolute difference per category|')
        )
    ggsave(plot=p, filename='binary_importance_by_category.pdf', 
                   limitsize=False, verbose = False)
    
    

    
def add_plot_helper_columns(df):
    assert ('annotation' in df.columns) == True, 'missing annotation column'
    df['anno_code'] = df.apply(lambda x: 
                                       BINARY_ANNOTATION_CODES[x['annotation']], axis=1)
    df['anno_category'] = df.apply(lambda x: 
                                       BINARY_ANNOTATION_CATEGORIES[x['annotation']], axis=1)
    df['anno_print'] = df.apply(lambda x: 
                                        ANNOTATION_CLEAN[x['annotation']], axis=1)
    
    return df 




'''
## for variance plot of binary annotations (relative importance) 
f = binary_annots.div(binary_annots.max(axis=1), axis=0)
f['pheno'] = f.index
ff = pd.melt(f, id_vars = ['pheno'], value_vars = list(set(f.columns) - set(['pheno'])))
ff = utils.add_plot_helper_columns(ff)
'''