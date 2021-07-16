#!/usr/bin/env python
#title       : optimize_specific_peaks.py
#description : Perform optimization for cell type specific peaks
#author      : Huamei Li
#date        : 27/08/2018
#type        : module
#version     : 2.7

#--------------------------------------------------------
# load python modules

from scipy import linalg

#--------------------------------------------------------
# load own modules

from modules.utils import *

#--------------------------------------------------------
# global setting

LOGS = log_infos() # logging informative

#--------------------------------------------------------

def extract_infos(profile, phenotypes, qvalues, fields, pi_score=1.0, exp_ratio=0.33):
    '''
    extract key infos for flow-up analysis
    :param profile: [pd.DataFrame] pure sample profile
    :param phenotypes: [pd.DataFrame] Phenotype classes file
    :param qvalues: [pd.DataFrame] q-value of cell type specific peaks
    :param fields: [list] target data field
    :param pi_score: [float] Pi-score value to descide cell type-specific peaks, default: 1.0
    :param exp_ratio: [float] peak expression concentration ratio, default: 0.33
    :return: infos [pd.DataFrame]
    
    '''
    (nrows, ncols), fold_changes, bool_lst = profile[fields].shape, [], []
    min_qvals = np.min(qvalues, axis=1)
    max_index = np.argmin(qvalues.values, axis=1)
    
    for idx, row in enumerate(profile[fields].values):
        index = max_index[idx]
        topv, others = row[index], np.delete(row, index)
        fchang = np.log2(topv / np.mean(others))
        fold_changes.append(fchang)
        bool_lst.append(topv * exp_ratio > np.max(others))
    
    infos = pd.DataFrame(
            np.array([-np.log10(min_qvals.values), max_index, fold_changes, bool_lst]).T,
            columns = ['Qvalue', 'TopIndex', 'FoldChange', 'Bool'],
            index = profile.index)
    
    profile = pd.concat([profile, infos], axis=1)
    profile = profile.loc[~(profile.Bool == 1)]
    profile['Score'] = profile.Qvalue * profile.FoldChange
    profile = profile[profile.Score >= pi_score ]
    LOGS.info('{} cell type specific peaks across {} cell types were identified'.format(profile.shape[0], phenotypes.shape[0]))

    return profile

def filter_peaks(profile, phenotypes, group_size):
    '''
    filter high confidence cell type specific peaks
    :param profile: [pd.DataFrame] pure sample profile and supplementary information, including q-value, top index, foldchange and specific score
    :param phenotypes: [pd.DataFrame] Phenotype classes file
    :param qvalues: [pd.DataFrame] q-value of cell type specific peaks
    :param group_size: [int] number of cell type specific peaks to consider for each phenotypes
    :return: sub_profile [pd.DataFrame]
    
    '''
    sub_profile, cellnames = None, phenotypes.index
    profile = profile.sort_values(by=['Score'], ascending=False)
     
    for idx, sample_idx in enumerate(profile.TopIndex.unique()):
        sub_tmp  = profile[(profile.TopIndex == sample_idx)]
        min_size = min(group_size, sub_tmp.shape[0])
        sub_tmp  = sub_tmp.iloc[0 : min_size]
        sub_profile = pd.concat([sub_profile, sub_tmp]) if idx else sub_tmp
    
    tmp_data  = sub_profile[phenotypes.index]
    tmp_data -= np.mean(tmp_data, axis=0)
    tmp_data  = tmp_data.divide(np.std(tmp_data, axis=0), axis=1)
    U, S, V   = linalg.svd(tmp_data, full_matrices=False)
    cond_val  = np.max(abs(S)) / np.min(abs(S))
    return sub_profile, cond_val

def optimize_peaks(profile, phenotypes, qvalues, min_group_size, max_group_size, pi_score=1.0, exp_ratio=0.33):
    '''
    optimize cell type peaks
    :param profile: [pd.DataFrame] pure sample profile
    :param phenotypes: [pd.DataFrame] Phenotype classes file
    :param qvalues: [pd.DataFrame] q-value of cell type specific peaks
    :param min_group_size: [int] minimum number of cell type specific peaks to consider for each phenotypes
    :param max_group_size: [int] maximum number of cell type specific peaks to consider for each phenotypes
    :param pi_score: [float] Pi-score value to descide cell type-specific peaks, default: 1.0
    :param exp_ratio: [float] peak expression concentration ratio, default: 0.33
    :return: signature_mat [pd.DataFrame]
    
    '''
    fields, cond_pre, size_pre = phenotypes.index, sys.float_info.max, 0
    profile = extract_infos(profile, phenotypes, qvalues, fields, pi_score, exp_ratio)
    group_size = range(min_group_size, max_group_size + 1)[::-1]
    del qvalues

    profile_bak = profile
    for idx, size in enumerate(group_size):
        sub_profile, cond_cur = filter_peaks(profile, phenotypes, size)
        if cond_cur > cond_pre: continue
        profile, cond_pre, size_pre = sub_profile, cond_cur, size
    
    sub_profile, cond_cur = filter_peaks(profile, phenotypes, size_pre)
    LOGS.info('Group size of each phenotype is {}, matrix condition number is {}'.format(size_pre, cond_cur))
    sub_profile = sub_profile[['chrom', 'start', 'end'] + phenotypes.index.tolist()]
    return sub_profile, profile_bak
