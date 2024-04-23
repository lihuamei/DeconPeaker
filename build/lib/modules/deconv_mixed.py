#!/usr/bin/env python
#title       : deconv_mixed.py
#description : Deconvolution of mixed samples based on pure profile. 
#author      : Huamei Li
#date        : 19/06/2018
#type        : module
#version     : 3.8

#-----------------------------------------------------
# load own modules

from modules.build_models import *
from modules.stat_plot    import *

#-----------------------------------------------------
# global setting

LOGS = log_infos() # logging informative

#-----------------------------------------------------

def normalize_profile(profile, method, fields, log=True, outfile=None):
    '''
    normalize pure profile to remove batch bias
    :param profile: [pd.DataFrame] profile matrix of peaks for all samples
    :param method: [str] specify normalize method
    :param fields: [list] data fields which need to be normalized
    :param log: [bool] log-transfer data or not, default: True
    :param outfile: [str] prefix name of output file, default: None
    :return: normalize profile [pd.DataFrame]
    
    '''
    norm_funcs = {
            'QN'   : quantile_norm      ,
            'DESeq': deseq_norm         ,
            'UQN'  : upper_quantile_norm,
            'PPM'  : ppm_norm           ,
            'TMM'  : tmm_norm
        }
    
    profile = norm_funcs[method](profile, fields)
    if log:
        profile[fields] = np.log2(1 + profile[fields])
    if outfile:
        profile.to_csv(outfile + '_{}_normalized.xls'.format(method), sep='\t', header=True, index=False)
    return profile

def filter_weakpeaks(profile):
    '''
    remove weak signal peaks and low variable peaks across different cell types
    :param profile: [pd.DataFrame] read counts of samples
    :return: pre-processed results
    
    '''
    colnames  = profile.columns[3 : ]
    lower_val = np.percentile(profile[colnames].values.flatten(), 50)
    profile   = profile.loc[~(np.sum(profile[colnames] < lower_val, axis=1) == len(colnames))]

    return profile

def cellspecificpeaks(profile, phenotype, kargs=None):
    '''
    obtained cell specific peaks of pure cells
    :param profile: [pd.DataFrame] pre-processed read counts
    :param phenotype: [pd.DataFrame] phenotype classes of all samples
    :param method: [str] select which method to find cell specific peaks
    :param thread: [int] thread number, default: None
    :param kargs: [int] other parameters, some of them may be used
    :return: 0
    
    '''
    sigmatrix, ctsp_peaks = find_marker_peaks(
            profile             , 
            phenotype           , 
            kargs.thread        , 
            kargs.score         , 
            kargs.min_group_size, 
            kargs.max_group_size,
            kargs.ratio         ,
            kargs.merge_replicates
        )
    fields = sigmatrix.columns[3 : ]
    if kargs.lib_strategy in ['RNA-Seq', 'Microarray']:
        bool_v, sigmatrix = True, sigmatrix[fields]
    else:
        bool_v = False
    bars(pd.DataFrame({
            'val': np.sum(ctsp_peaks[phenotype.index], axis=0), 
            'lab': phenotype.index.tolist()
        }), os.path.join(kargs.outdir, kargs.prefix + '_cstps_counts'), platform=kargs.lib_strategy)

    outfile = os.path.join(kargs.outdir, kargs.prefix + '_signature_matrix.xls')
    sigmatrix.to_csv(outfile, sep='\t', header=True, index=bool_v)
    
    outfig = os.path.join(kargs.outdir, kargs.prefix + '_signature_heatmap')
    if sigmatrix.shape[0] <= 10000: cluster_heatmap(sigmatrix[fields], outfig)
    return sigmatrix

def deconvcells(mixsamples, sigprofile, lib_strategy=None, pvalue=False, method='SIMPLS', norm=None, outdir='./'):
    '''
    mixed samples will be deconvolved based on the signal of the cell-specific peaks of the pure cells
    :param mixsamples: [pd.DataFrame] multiple-mixed-samples singal matrix file
    :param sigprofile: [pd.DataFrame] pure cell profile file
    :param lib_strategy: [str] a string indicating the type of the profile measurements, default: None
    :param pvalue: [bool] estimate p-value or not, default: False
    :param method: [str] deconvolution method, including SIMPLS and RSIMPLS, default: SIMPLS
    :param norm: [str] normalize method, default: None
    :param outdir: [str/dir] output directory, default: ./
    :return: deconvoluted results file
    
    '''  
    mixsamples, sigprofile = intersect(mixsamples, sigprofile, lib_strategy)
    mixsamples, sigprofile = mixsamples[mixsamples.columns[3 : ]], sigprofile[sigprofile.columns[3 : ]]
    if norm: 
        mixsamples = normalize_profile(
                mixsamples        , 
                norm              , 
                mixsamples.columns, 
                log = False       , 
                outfile = None
            )       
    deconv_results = deconv(mixsamples, sigprofile, method=method, pvalue=pvalue)
    outfile = os.path.join(outdir, 'deconPeaker-Results.xls')
    deconv_results.to_csv(outfile, sep='\t', index=True, header=True)
    
    outfig = os.path.join(outdir, 'deconPeaker-Results')
    stack_bars(deconv_results[sigprofile.columns], outfig)
    return deconv_results
