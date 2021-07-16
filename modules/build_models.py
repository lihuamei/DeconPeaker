#!/usr/bin/env python
#title       : build_models.py
#description : Build models for Deconvoluting mixed samples based on pure profile. 
#author      : Huamei Li
#date        : 19/06/2018
#type        : module
#version     : 3.8

#-----------------------------------------------------
# load python modules

from rpy2.robjects   import r, pandas2ri

#-----------------------------------------------------
# load own modules

from modules.find_markers import *

#-----------------------------------------------------
# global setting

pandas2ri.activate()
LOGS = log_infos() # logging informative

#--------------------------------------------------------------------------
# SIMPLS deconvolution method

def SIMPLS(Y, X, method, pvalue=False):
    '''
    :param Y: [pd.DataFrame] mixed sample profile
    :param X: [pd.DataFrame] pure cells profile
    :param method: [str] deconvolution method, include PCR and SIMPLS
    :return: deconvolution results
    
    '''
    r.assign('Y', Y); r.assign('X', X); r.assign('method', method); r.assign('pvalue', pvalue)
    
    r('''
        source('modules/simpls_deconv.r')
        results <- simpls_deconv(Y, X, method, pvalue)
    ''')
    rmse, coeff, r2, pvalue = r('results$rmse')[0], r('results$coeffs'), r('results$R2')[0], r('results$pval')[0]
    infos  = [r2, rmse, pvalue]
    return np.append(coeff, infos)

#--------------------------------------------------------------------------
# RSIMPLS deconvolution method

def RSIMPLS(Y, X):
    '''
    :param Y: [pd.DataFrame] mixed sample profile
    :param X: [pd.DataFrame] pure cells profile
    :return: deconvolution results

    '''
    cur_path = os.getcwd()
    os.chdir('modules')
    eng = matlab_engine().start_matlab()
    
    Y_fp, X_fp = create_tmp_files(['Y_tmp', 'X_tmp'])
    Y.to_csv(Y_fp, sep='\t', header=True, index=True)
    X.to_csv(X_fp, sep='\t', header=True, index=True)
    Y_fp.close(); X_fp.close()
    
    results = eng.rsimpls_deconv(Y_fp.name, X_fp.name, nargout=1)
    os.chdir(cur_path)
    return np.array(results)

#--------------------------------------------------------------------------
# SIMPLS: An alternative approach to partial least squares regression

def deconv(Y, X, method='SIMPLS', pvalue=False):
    '''
    deconvolution using a SIMPLS strategy
    :param Y: [pd.DataFrame] mixed sample profile
    :param X: [pd.DataFrame] pure cells profile
    :param method: [str] deconvolution method, including SIMPLS and RSIMPLS, default: SIMPLS
    :param pvalue: [bool] estimate P-value or not, default: False
    :return: coeffs [np.array]
    
    '''
    mixsnames, purecells, deconv_results = Y.columns, X.columns, [] 
    
    if method == 'RSIMPLS':
        deconv_results = RSIMPLS(Y, X)
    else:
        for sn in mixsnames: deconv_results.append(SIMPLS(Y[sn], X, method, pvalue))
        #deconv_results = [ SIMPLS(Y[sn], X, method, pvalue) for sn in mixsnames ]
    deconv_results = pd.DataFrame(
            deconv_results,
            columns = np.append(purecells, ['Rsquared', 'RMSEP', 'P.value']),
            index = mixsnames
        )
    return deconv_results
