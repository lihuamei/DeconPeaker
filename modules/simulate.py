#!/usr/bin/env python
#title      : simulate.py
#decription : Simulate reads with fixed proportion
#author     : Huamei Li
#date       : 21/06/2018
#type       : module
#version    : 3.8

#-----------------------------------------------------
# load own modules

from modules.bams import *
from modules.random_proportion import *

#-----------------------------------------------------
# global setting
    
LOGS = log_infos() # logging informative

#-----------------------------------------------------
def sim_from_bam(pure_bams, readcnts, kargs):
    '''
    sampling reads from BAMs of pure cells
    :param pure_bams: [str/path] bam of pure cells
    :param readcnts: [int] total number of reads after sampling 
    :param kargs: [dict] total input parameters of simulation step
    :return: tmpbams [list/files]
    
    '''
    tmpbams = []
    for idx, bam in enumerate(pure_bams):
        bam   = bam[0] if isinstance(bam, np.ndarray) else bam
        infos = kargs.pure_infos[kargs.pure_infos['data'] == bam]
        cell  = infos['cellname'][0]
        prob  = readcnts[cell] / float(kargs.cellcounts[cell])
        tmpfp = create_tmp_files('.bam')[0].name
        cmds  = 'samtools view -s {0} -b {1} -@ 4 > {2}'.format(prob, bam, tmpfp)
        syscmd_run(cmds)
        tmpbams.append(tmpfp)
    return tmpbams

def adjust_samplecounts(readcounts, proportion, purecells, kargs):
    '''
    adjust read counts for sampling if simulated read counts exceed total read counts, 
    :param readcounts: [int] sampling reade counts which defined by users
    :param proportion: [list] proportion of each cell
    :param purecells: [list] pure cells needed for mixing
    :param kargs: [dict] all of input parameters
    :return: adj_readcnts [list] read counts of each cell for sampling
    
    '''
    adj_readcnts = [ int(readcounts * ratio) for ratio in proportion ]
    totalcnts = [ kargs.cellcounts[cell] for cell in purecells ]
    while 1:
        flag = 0
        for idx, rawcnt in enumerate(totalcnts):
            if adj_readcnts[idx] > rawcnt: 
                flag += 1
                break
        if not flag: break
        readcounts *= 0.99
        adj_readcnts = [ int(readcounts * ratio) for ratio in proportion ]
    LOGS.warn('Sub-sampling read counts for {} with proportions {}: {}'.\
            format(', '.join(purecells), ', '.join(map(str, proportion)), ', '.join(map(str, adj_readcnts))))
    return { cell : cnts for cell, cnts in zip(purecells, adj_readcnts) }

def simulate_mixture(kargs):
    '''
    Generate mixture samples with simulated pure cell proportions
    :param kargs: [dict] all of input parameters
    :return: 0
    
    '''
    readcounts = kargs.readcounts
    for idx, rowinfos in kargs.mixture.iterrows():
        proportion, purecells = rowinfos[rowinfos > 0], rowinfos[rowinfos > 0].index
        adj_readcnts = adjust_samplecounts(readcounts, proportion, purecells, kargs)
        datafiles    = kargs.pure_infos[kargs.pure_infos['cellname'].isin(purecells)]['data'].values
        tmpfiles = multi_process(
                    datafiles              , 
                    sim_from_bam           ,
                    kargs.thread           , 
                    readcnts = adj_readcnts,
                    kargs = kargs
                ) # simulation reads of each cell with fix proportion
        prefix = kargs.prefix + '_' + rowinfos.name
        outdir = os.path.join(kargs.outdir, str(len(adj_readcnts)), '_'.join(adj_readcnts.keys()))		
        merge_bams(tmpfiles, prefix, outdir=outdir, threads=kargs.thread)
        [ os.remove(fil) for fil in tmpfiles ]
    return 0

def get_pure_datafiles(kargs):
    '''
    get pure data files for samping
    :param kargs: [dict] all of input parameters
    :return: sub_files [list] 
    
    '''
    maxidx = np.argmax(np.sum(kargs.mixture > 0, axis=1))
    tmpidx = np.where(kargs.mixture.loc[maxidx, :] > 0)
    purecells = kargs.mixture.columns[tmpidx]
    purecells = list(set(purecells))
    datafiles = dict(zip(kargs.pure_infos.cellname, kargs.pure_infos.data))
    sub_files = [ datafiles[cell] for cell in purecells ]
    return sub_files, purecells

def bed2bam(kargs):
    '''
    convert BED file to BAM file
    :param kargs: [dict] all of input parameters
    :return: kargs
    
    '''
    tmpbamfiles = []
    for bed in kargs.pure_infos.data:
        tmpfp = create_tmp_files('.bam')[0].name
        cmds  = 'bedToBam -i {} -g {} > {}'.format(bed, kargs.genome, tmpfp)
        syscmd_run(cmds)
        tmpbamfiles.append(tmpfp)
    kargs.pure_infos.data = tmpbamfiles
    return kargs

def multi_simulate_bams(kargs):
    '''
    multi-processes simulate mixed cell samples fromm BAM files of pure cells
    :param kargs: [dict] all of input parameters
    :return: 0
    
    '''
    if kargs.format == 'BED': kargs = bed2bam(kargs)
    bamfiles, purecells = get_pure_datafiles(kargs)
    LOGS.info('Calculating the number of original reads of pure cell samples')
    readcounts = multi_process(
            bamfiles      , 
            get_readcounts, 
            kargs.thread
        )
    kargs.cellcounts = dict(zip(purecells, readcounts))
    LOGS.info('Proportional sampling of tags......')
    kargs.mixture = kargs.mixture[purecells]
    simulate_mixture(kargs)
    return 0
