#!/usr/bin/env python
#title       : norm_counts.py 
#description : Normalize read counts of cell types based on chromatin accessibility. 
#author      : Huamei Li
#date        : 29/05/2018
#type        : module script
#version     : 3.8

#-----------------------------------------------------
# load own python modules

from modules.utils  import *

#-----------------------------------------------------
# global setting

LOGS = log_infos() # logging informative

#-----------------------------------------------------

def quantile_norm(profile, fields):
    '''
    normalize profile by quantile method
    :param profile: [pd.DataFrame] pure sample profile, row=peaks, column=cell
    :param fields: [list] field names of profile which need to be normalize
    :return: nomalized profile
    
    '''
    X = profile[fields]
    quantiles = np.mean(np.sort(X, axis=0), axis=1)
    ranks = np.apply_along_axis(stats.rankdata, 0, X)
    rank_indices = ranks.astype(int) - 1
    profile[fields] = pd.DataFrame(quantiles[rank_indices], index=profile.index)
    return profile

def deseq_norm(profile, fields):
    '''
    normalize profile by deseq method
    :param profile: [pd.DataFrame] pure sample profile, row=peaks, column=cell
    :param fields: [list] field names of profile which need to be normalize
    :return: nomalized profile
    
    '''
    pse_ref = np.mean(np.log(profile[fields] + 1), axis=1)
    ratios  = np.log(profile[fields] + 1).sub(pse_ref, axis=0)
    factors = np.exp(np.median(ratios, axis=0))
    profile[fields] = profile[fields] / factors
    return profile

def upper_quantile_norm(profile, fields, q=0.75):
    '''
    normalize profile by upper quantile method
    :param profile: [pd.DataFrame] pure sample profile, row=peaks, column=cell
    :param fields: [list] field names of profile which need to be normalize
    :param q: [float] quantile, default: 0.75
    :return: nomalized profile
    
    '''
    lib_size = profile[fields].sum()
    tmp_val  = profile[fields].T.div(lib_size, axis=0).T
    tmp_val  = tmp_val.quantile(q)
    factors  = tmp_val.div(np.exp(np.mean(np.log(tmp_val))))

    profile[fields] = profile[fields].div(factors, axis=1)
    return profile

def ppm_norm(profile, fields):
    '''
    normalize read counts by PPM method
    :param readcounts: [pd.DataFrame] read counts data frame, row=peaks, column=cell
    :param fields: [list] field names of readcounts which need to be normalize
    :return: nomalized readcounts
   
    '''
    profile[fields] = profile[fields] / np.sum(profile[fields]) * 1e6
    return profile

def tmm_norm(profile, fields):
    '''
    normalize read counts by TMM method
    :param readcounts: [pd.DataFrame] read counts data frame, row=peaks, column=cell
    :param fields: [list] field names of readcounts which need to be normalize
    :return: nomalized readcounts

    '''
    r.assign('data', profile[fields])
    r('''
        library(edgeR)
        factors <- calcNormFactors(as.matrix(data), method = 'TMM')
    ''')
    factors = r('factors')
    profile[fields] = profile[fields].div(factors, axis=1)
    return profile

