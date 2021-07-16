#!/usr/bin/env python
#title       : deconPeaker.py 
#description : Identification and deconvolution of cell types based on chromatin accessibility. 
#author      : Huamei Li
#date        : 29/05/2018
#type        : main script
#version     : 3.8

#-----------------------------------------------------
# load python modules

import numexpr
from time import time

#-----------------------------------------------------
# load own python modules

from modules.deconv_mixed import *
from modules.simulate     import *
from modules.peaks        import *
from modules.parse_opts   import parse_opts

#-----------------------------------------------------
# global setting

LOGS = log_infos() # logging informative

numexpr.set_num_threads(numexpr.detect_number_of_cores())

#-----------------------------------------------------

def preprocess():
    '''
    Create chromatin accessibility profile for pure samples, 
    this step is the basis for subsequent specific cell type identification and mixture deconvolution.
    :return: 0

    '''
    LOGS.info('Loading peaks......')
    merged_files = mp_read_peaks(ARGS.infos.PEAK.values.tolist(), kargs=ARGS)
    
    LOGS.info('Chosing a list of non-overlapping, maximally significant peaks')
    nonovp_peakfil = remove_redundant_peakfile(
            merged_files              , 
            ARGS.infos.CELL.unique()  , 
            ARGS.prefix               , 
            ARGS.outdir
        )
    
    if ARGS.offset:
        LOGS.info('Filtering out the peaks nearby the TSS (+/-{} bps)'.format(ARGS.offset))
    nonovp_peakfil = remove_peaks_nearbytss(nonovp_peakfil, ARGS.hg_genome, offset=ARGS.offset)
    
    peaknum = get_line_number(nonovp_peakfil)
    LOGS.info('Counting the number of fragments from each sample falling into each of {} peaks'.format(peaknum))
    multi_get_reads(
            ARGS.infos.BAM.values.tolist(), 
            nonovp_peakfil                , 
            ARGS                          , 
            prefix = ARGS.prefix + '_reference_count_matrix',
            outdir = ARGS.outdir          ,
            bg = False
        )
    phenotype = os.path.join(ARGS.outdir, ARGS.prefix + '_phenotype_classes.xls')
    LOGS.info('Writing phenootype class into {}'.format(phenotype))
    write_phenotype_file(ARGS.infos, phenotype)
    LOGS.info('Preprocess finished')
    return 0

def findctsps():
    '''
    find cell type specific peaks for pure samples
    :return: 0

    '''
    profile    = load_profile(ARGS.profile, ARGS.lib_strategy)
    phenotypes = load_phenotypes(ARGS.phenotype)
    LOGS.info('Loaded {} peaks and {} samples'.format(profile.shape[0], profile.shape[1] - 3))

    ARGS.prefix = os.path.basename(ARGS.profile).rsplit('.', 1)[0]
    LOGS.info('Normalizing pure profile by {} method to remove batch effects'.format(ARGS.norm))
    profile_norm = normalize_profile(
            profile              , 
            ARGS.norm            , 
            profile.columns[3 : ],
            False                
        )
    
    profile_norm = filter_weakpeaks(profile_norm)
    curcnts, diffcnts = profile_norm.shape[0], profile.shape[0] - profile_norm.shape[0]
    LOGS.info('Filtering out {} peaks and {} peaks have been remained'.format(diffcnts, curcnts))
    del profile

    LOGS.info('Identifying cell specific peaks accross pure cell profile')
    markerpeaks = cellspecificpeaks(profile_norm, phenotypes, ARGS)

    LOGS.info('Final number of cell specific peaks is {}'.format(markerpeaks.shape[0]))
    return 0

def deconvolution():
    '''
    Based on pure cell type information, a deconvolution strategy was used to 
    calculate the proportion of possible cell types in the mixed sample.
    
    '''
    if ARGS.format == 'BAM':
        LOGS.info('Counting the number of reads/fragments in each peak in the bam files by featureCounts')
        ARGS.mixture = multi_get_reads(
            ARGS.infos.BAM.values.tolist() ,
            ARGS.pure                      ,
            ARGS                           ,
            prefix = 'mixed_sample_profile',
            outdir = ARGS.outdir           ,
            bg = False
        )
    mixsamples, sigprofile = load_profile([ARGS.mixture, ARGS.pure], ARGS.lib_strategy)
    LOGS.info('Deconvolving......')
    results = deconvcells(
            mixsamples       , 
            sigprofile       , 
            ARGS.lib_strategy,
            ARGS.pvalue      ,
            ARGS.method      ,
            ARGS.norm        ,
            ARGS.outdir      
        )
    LOGS.info('Showing deconPeaker results: ')
    print(results)
    return 0

def simulate():
    '''
    simulate mixture cell sampels by sampling reads of pure cells
    
    '''
    LOGS.info('Simulating......')
    ARGS.mixture = random_proportions(
            ARGS.pure_infos, 
            ARGS.mixture   , 
            ARGS.rep_counts, 
            ARGS.prefix    , 
            ARGS.outdir
        )
    multi_simulate_bams(ARGS) # simulate data with specified proportions
    LOGS.info('Okay!')
    return 0

def run():
    '''
    deconPeaker main funcion and contains four sections, pureprofile, identifycells, deconvolution and simulation, respectivly
    :return: stat [int] 

    '''
    global ARGS
    ARGS, start  = parse_opts(), time()
    funcs, modes = [preprocess, findctsps, deconvolution, simulate], \
                   ['preprocess', 'findctsps', 'deconvolution', 'simulation']
    ARGS.tmpdir = tmpdir = __import__('tempfile').mkdtemp('_deconPeaker')
    try:    
        stat = funcs[modes.index(ARGS.sub_parser)]()
    finally:
        __import__('shutil').rmtree(tmpdir) if os.path.exists(tmpdir) else 0
    LOGS.info('Elapsed time is {} seconds'.format(time() - start))
    return stat

if __name__ == '__main__':
    sys.exit(run())

