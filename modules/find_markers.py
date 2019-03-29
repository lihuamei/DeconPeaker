#!/usr/bin/env python
#title       : find_marker.py
#description : Filtered cell specific peaks from different cells
#author      : Huamei Li
#date        : 27/06/2018
#type        : module
#version     : 2.7

#--------------------------------------------------------
# load own modules

from normalize_methods       import *
from optimize_specific_peaks import *
from lm_reg import get_cell_specific_pvals

#--------------------------------------------------------
# global setting

LOGS = log_infos() # logging informative

#--------------------------------------------------------

def find_marker_peaks(profile, phenotypes, *kargs):
    '''
    find marker cell type specific peaks accross cells
    :param profile: [pd.DataFrame] pure sample profile
    :param phenotypes: [pd.DataFrame] replicate file which contains cell positon in the dataframe
    :return: sigmatrix [pd.DataFrame]
    
    '''
    fields  = profile.columns[3 : ]
    threads, pi_score_cutoff, min_group_size, max_group_size, exp_ratio, merge_method = kargs
    
    qvalues = get_cell_specific_pvals(
            profile   , 
            phenotypes, 
            fields    , 
            threads = threads
        )
    
    pdb.set_trace()
    profile = profile[~profile.index.duplicated(keep='first')]
    qvalues = qvalues[~qvalues.index.duplicated(keep='first')]
    profile = merged_replicates(profile, phenotypes, method=merge_method).loc[qvalues.index]
    
    LOGS.info('{} cell type specific peaks across {} cell types were identified'.format(profile.shape[0], phenotypes.shape[0]))
    
    LOGS.info('Performing optimization for cell type specific peaks')
    sigmatrix = optimize_peaks(
            profile        , 
            phenotypes     , 
            qvalues        , 
            min_group_size , 
            max_group_size ,
            pi_score_cutoff,
            exp_ratio      
        )
    return sigmatrix
