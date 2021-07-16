#!/usr/bin/env python
#title       : lm_reg.py 
#description : Get pvalues of peaks of each cell type
#author      : Huamei Li
#date        : 29/05/2018
#type        : main script
#version     : 3.8

#-----------------------------------------------------
# load own modules

from modules.utils import *

#-----------------------------------------------------
# global setting

LOGS = log_infos() # logging informative

#-----------------------------------------------------
def design_bin(phenotypes):
    '''
    create design matrix, x(i, k) = 1 if sample belongs to cell type k and x(i, k) = 0 otherwise.
    :param phenotypes: [pd.DataFrame] phenotype informative
    :return: design_matrix [pd.DataFrame] and the number of samples [np.array]
    
    '''
    design_matrix = __import__('copy').deepcopy(phenotypes)
    design_matrix[ design_matrix == 2 ] = 0
    return design_matrix

def contranst_mat(ncells):
    '''
    get all the contrasts that is needed to calculate the p-value
    :param ncells: [int] cell counts
    :return: X [np.matrix]
    
    '''
    nrows, ncols = (ncells - 1) * ncells, ncells
    X, flag, celltype = np.zeros((nrows, ncols)), 0, []
    
    for i in range(nrows):
        X[i, int(i / (ncols - 1))] = 1
        if flag == i / (ncols - 1): flag += 1
        X[i, flag] = -1
        if (flag + 1) / ncells:
            flag = 0
        else:
            flag += 1
        celltype.append(i / (ncols - 1))
    
    return {'mat': X, 'type': celltype }

def multi_lmreg(datalst, X, XX_inv_X, XX_inv, contrasts_dict):
    '''
    build multiple regression and calculate pvalues
    :param datalst: [pd.DataFrame in list] dependent variable
    :param X: [pd.DataFrame] design matrix, indenpendent variables
    :param XX_inv: [pd.DataFrame] XX_inv = inv((X'X))X'
    :param contrasts_dict: [np.matrix] contrast dictionary, include contrast matrix, pairwise sum of sample counts
    :return pv_infos [list in list]
   
    '''
    pvalst, (nrows, celltypecnt) = [], X.shape
    data, cellnames, peaknames = datalst[0].values, X.columns.tolist(), datalst[0].index
    obs_cnts = np.sum(X, axis=0)
    contrast, cell_cluster = contrasts_dict['mat'], np.array(contrasts_dict['type'])
    
    for idx, Y in enumerate(data):
        betas = np.dot(XX_inv_X, Y)
        resid = Y - np.dot(X, betas)
        sigma = np.dot(resid ** 2, X) / (obs_cnts)

        beta_differ = np.dot(contrast, betas)
        var_betahat = np.dot(np.diag(sigma * (obs_cnts - 1)), XX_inv)
        beta_se_val = np.sqrt(np.diagonal(np.dot(np.dot(contrast, var_betahat), contrast.T)) / (contrasts_dict['cnts'] - 2))
    
        pvals = stats.t.sf(beta_differ / beta_se_val, np.dot(np.abs(contrast), obs_cnts) - 2)
        pvals = [ np.max(pvals[cell_cluster == i]) for i in range(celltypecnt) ]
        pvalst.append(pvals)
    return pvalst

def get_cell_specific_pvals(profile, phenotypes, fields, threads=1):
    '''
    get cell type specific peak pvalues
    :param profile: [pd.DataFrame] pure sample profile matrix
    :param phenotypes: [pd.DataFrame] phenotype informative
    :param fields: [list] profile fileds for buld regression model
    :param threads: [int] threads number, default: 1
    :return: qvalues [pd.DataFrame]
    
    '''
    design_matrix = design_bin(phenotypes).T
    peak_profiles = np.array_split(profile[fields], threads, axis=0)
    XX_inv   = np.linalg.inv(np.dot(design_matrix.T, design_matrix))
    XX_inv_X = np.dot(XX_inv, design_matrix.T)
    contrasts_dict = contranst_mat(design_matrix.shape[1])
    contrasts_dict['cnts'] = np.dot(np.abs(contrasts_dict['mat']), np.sum(design_matrix, axis=0))
    
    pval_infos = multi_process(
            peak_profiles       , 
            multi_lmreg         , 
            threads             , 
            X = design_matrix   ,
            XX_inv_X = XX_inv_X ,
            XX_inv = XX_inv     ,
            contrasts_dict = contrasts_dict
        )
    pval_infos = pd.DataFrame(pval_infos, columns=design_matrix.columns)
    pval_infos = pval_infos.fillna(1)
    r.assign('pvalues', pval_infos)
    del pval_infos
    
    r('''
        pvalues = as.matrix(pvalues)
        qvalues = matrix(p.adjust(as.vector(pvalues), method='BH'), ncol=ncol(pvalues))                    
    ''')
    qvalues = pd.DataFrame(r('qvalues'), columns=phenotypes.index, index=profile.index)
    return qvalues
