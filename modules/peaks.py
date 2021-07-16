#!/usr/bin/env python
#title      : filter_peaks.py
#decription : Take various treatment on the peaks
#author     : Huamei Li
#date       : 02/06/2018
#type       : module
#version    : 3.8

#-----------------------------------------------------
# load python modules

from scipy import io

#-----------------------------------------------------
# load own modules

from modules.utils import *

#-----------------------------------------------------
# global setting

LOGS = log_infos()

if sys.platform == 'win32':
	pass
else:
	from bx.intervals.intersection import Intersecter, Interval

#----------------------------------------------------
def select_peaks(peakfil, merged_fp, cut_qval=2, cellname=None):
    '''
    filter out weak signal peaks, and then adjust to the uniform width if uniform=True
    :param peakfil  : [str/file] enriched regions called by MACS2
    :param merged_fp: [hanlde] file handler
    :param cut_pval : [float] critical qvalue between weak and stringent peaks, defaullt: qvalue = 0.01 (-log10(q)=2)
    :param cellname : [str] cell name, default: None
    :return: merged_fp [handler] merged file handler and peaks [pd.DataFrame]
      
    '''
    peaks = pd.read_csv(peakfil, sep='\t', names=NARROWS_NAMES)
    peaks = peaks[peaks['qValue'] >= cut_qval]
    peaks = peaks[peaks['chrom'].isin(CHROMS)]
    
    peaks['cellname'], peaks['width'] = cellname, peaks['end'] - peaks['start']
    summits_pos = peaks['start'] + peaks['summit2PeakDist']
    peaks['start'], peaks['end'] = summits_pos - 250, summits_pos + 250
    
    peaks = peaks.sort_values(by=['chrom', 'start', 'end'], ascending=True)
    peaks.to_csv(merged_fp, mode='a', sep='\t', header=False, index=False)
    return merged_fp, peaks

def single_read_peaks(tasks, kargs):
    '''
    single-process filters peaks according to preset conditions
    :param tasks: [list] the amount of tasks that a single-process needs to perform
    :param kargs: [dict] other param infos, including cut_pval, cut_qval and blacklist (blacklist file path)
    :rerurn: [str] temp_path
    
    '''
    merged_file = create_tmp_files(mode='a')[0]
    for tk in tasks:
        peak_fil = filter_blacklist(tk, kargs.blacklist)
        cellname = kargs.infos[kargs.infos['PEAK'] == tk]['CELL'].values.tolist()[0]
        merged_file = select_peaks(
                peak_fil               ,
                merged_file            ,
                cut_qval = kargs.qvalue,
                cellname = cellname    ,
            )[0] # select stringent peaks
    merged_file.close()
    return [ merged_file.name ]

def filter_blacklist(peak_file, blacklist):
    '''
    filter out peaks which in the blacklist
    :param peak_file: [str/file] merged temporary file
    :param blacklist: [str] blacklist file which download from ENCODE
    :return: non-blacklist peak file
    
    '''
    peak_file = decomp_gzfile(peak_file)[0]
    noblacklist = create_tmp_files('.no_blacklist')[0].name
    
    cmd = 'sub_tools/bedtools intersect -v -a {} -b {} > {}'.format(
            peak_file ,
            blacklist ,
            noblacklist
        )
    syscmd_run(cmd, rm=peak_file) # remove peaks which in the blacklist
    return noblacklist

def find_ovp_lsts(container, idx=6):
    '''
    find overlap for two list which have been sorted
    :param container: [list, list] peak region information
    :param idx: [int] peak score position, default: 6
    :return: container
    
    '''
    lst1, lst2 = container
    if lst1[0] == lst2[0] and lst2[1] < lst1[2]:
        flag = 0 if lst2[idx] > lst1[idx] else 1
        container.pop(flag)
    return container

def remove_redundant_sorted_peaks(peaks, fp):
    '''
    quickly remove redundant peaks for pandas dataframe
    :param peaks: [pd.DataFrame] peaks informative
    :param fp: [handler] file handler which peaks need to be save
    :return: non-redundant peak file handler
              
    '''
    columns, peaks = peaks.columns, peaks.values
    (nrows, ncols), container = np.shape(peaks), []
    
    for idx, peak in enumerate(peaks):
        container.append(peak)
        if len(container) / 2 != 1: continue
        container = find_ovp_lsts(container) 

        if len(container) < 2: continue
        fp.write('\t'.join(map(str, container[0])) + '\n')
        container.pop(0)
    
    fp.write('\t'.join(map(str, container[0])) + '\n')
    return fp

def remove_redundant_peakfile(peakfil, cells, prefix, outdir):
    '''
    choose a list of non-overlapping peaks by specified strategy
    :param peakfil: [str/file] path of the peak file
    :param cells: [list] unique cell names of pure samples
    :param prefix: [str] prefix name of non-redundant peaks output file
    :param outdir: [str/dir] output directory
    :return: path of non-overlapping peak file
    
    '''
    tmp_overlap = create_tmp_files('.merged_nonoverlap_peaks.bed')[0]
    cell_map, container = { cell : idx for idx, cell in enumerate(cells) }, []
    with open(peakfil, 'rb') as fp_r:
        for idx, line in enumerate(fp_r):
            line_sep = line.decode().strip().split('\t')
            (chrom, start, end), score, cell = line_sep[0 : 3], line_sep[6], line_sep[10]
            start, end, score = int(start), int(end), float(score)
            container.append([chrom, start, end, score])

            if len(container) / 2 != 1: continue
            container = find_ovp_lsts(container, idx=3)
            if len(container) < 2: continue        
            tmp_overlap.write('\t'.join(map(str, container[0][0 : -1])) + '\n')
            container.pop(0)
    
    tmp_overlap.write('\t'.join(map(str, container[0][0 : -1])) + '\n')
    tmp_overlap.close()
    return tmp_overlap.name

def find_overlap_dataframes(query, hits):
    '''
    find overlap between sorted query and hits regions
    :param query: [pd.DataFrame] query peaks
    :param hits: [pd.DataFrame] hits regions
    :return: flags [array]
     
    '''
    query_idx, hits_labs = [], []
    tree_hash = { chrom : Intersecter() for chrom in hits.chrom.unique() }
    null_res  = [ tree_hash[chrom].add_interval(Interval(start, end)) for chrom, start, end in hits.values ]
    
    for idx, line in enumerate(query.values):
        chrom, start, end = line[0 : 3]
        overlaps = tree_hash[chrom].find(start, end)
        if not overlaps: continue
        for ovp in overlaps:
            hits_labs.append(chrom + ':' + str(ovp.start) + '-' + str(ovp.end))
            query_idx.append(idx)
    
    tmp_labels = hits.chrom + ':' + hits.start.astype(str) + '-' +  hits.end.astype(str)
    hit_indexs = tmp_labels[tmp_labels.isin(hits_labs)].index.values
    return np.array(query_idx), hit_indexs

def remove_peaks_nearbytss(peakfile, hgfile, offset, downstream=500):
    '''
    remove peaks which nearby TSS 
    :param peakfile: [str/file] peak file, which contains chromosomes start and end
    :param hgfile: [str/file] TSS postion file
    :param offset: [int] peaks falling into [TSS - offset, TSS + offset] need to be removed
    :param downstream: [int] offset distance of downstream, default: 500
    :return: non-redundant file
    
    '''
    columns  = ['chrom', 'start', 'end']
    tssinfos = pd.read_csv(hgfile, sep='\t', header=0)
    
    tssinfos.loc[tssinfos.strand == '+', 'start'] -= offset
    tssinfos.loc[tssinfos.strand == '+', 'end'] += downstream
    tssinfos.loc[tssinfos.strand == '-', 'start'] -= downstream
    tssinfos.loc[tssinfos.strand == '-', 'end'] += offset

    tssinfos.start, tssinfos.end = tssinfos.start - offset, tssinfos.end + offset
    tssinfos.loc[tssinfos['start'] < 0, 'start'] = 0

    peaks = pd.read_csv(peakfile, sep='\t', names=columns)
    flagindex = find_overlap_dataframes(peaks, tssinfos[columns])[0]
    peaks = peaks.drop(flagindex).reset_index(drop=True)
    peaks.to_csv(peakfile, sep='\t', header=False, index=False)
    return peakfile

def mp_read_peaks(peak_files, kargs):
    '''
    multi-processes to filter peaks, and then write into temporary files
    :param peak_files: [list] the corresponding peak file for each sample
    :param kargs: [dict] other param infos, including cut_pval, cut_qval , nth (thread number) and blacklist (blacklist file path)
    :return: merged_files [str/file] sorted peak file
     
    '''
    peak_files = list(set(peak_files)) if isinstance(peak_files, list) else [ peak_files ]
    merged_files = multi_process(
            peak_files       , 
            single_read_peaks, 
            kargs.thread     , 
            kargs = kargs
        )
    [ catfiles(fil, merged_files[0], rm=fil) for index, fil in enumerate(merged_files) if index > 0 ]
    merged_sorted = create_tmp_files('.merged_sorted_peaks.bed')[0].name
    sortfiles(merged_files[0], merged_sorted, rm=merged_files[0]) # quick sort peak file
    merged_files = merged_sorted
    return merged_files

