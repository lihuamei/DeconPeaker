#!/usr/bin/env python
#title      : bams.py
#decription : Read and parse BAM format file
#author     : Huamei Li
#date       : 11/06/2018
#type       : module
#version    : 3.8

#-----------------------------------------------------
# load python modules

import subprocess

#-----------------------------------------------------
# load own modules

from modules.utils import *

#-----------------------------------------------------
# global setting

if sys.platform == 'win32':
	pass
else:
	import pysam

LOGS = log_infos()

#-----------------------------------------------------

def convert_saf(bedfil, fg=True, upstream=500000, downstream=500000):
    '''
    convert bed to SAF format
    :param bedfil: [str/file] bed format file
    :param fg: [bool] compute foreground readcounts of cells or background, default: True
    :param upstream: [int] upstream of the peaks, default: 500000
    :param downstream: [int] downstream of the peaks, default: 500000
    :return: SAF file path
    
    '''
    suffix = '_foreground.saf' if fg else '_background.saf'
    safile = create_tmp_files(os.path.basename(bedfil) + suffix, fixnames=True)[0]
    with open(bedfil, 'rb') as fp_r:
        for line in fp_r:
            chrom, start, end = line.decode().strip().split('\t')[0 : 3]
            if chrom == 'chrom': continue
            start, end = map(int, [start, end])
            if not fg:
                start = 1 if start - upstream < 0 else start - upstream
                end  += downstream
            start, end = map(str, [start, end])
            label = chrom + '.' + start + '.' + end
            safile.write('\t'.join([label, chrom, start, end, '.']) + '\n')
    safile.close()
    return safile.name

def multi_get_reads(bamfils, bedfil, kargs, prefix='pure_cells_readcounts', outdir=None, bg=True):
    '''
    multi-process to get read counts for candidate peaks
    :param bamfils: [list] bam files
    :param bedfil: [str/file] non-redundant peak list file
    :param kargs: [dict] input parameters by users
    :param prefix: [str] prefix name of output read counts file, default: pure_cells_readcounts
    :param outdir: [str/dir] output directory, default: None
    :param bg: [bool] extract background read counts or not, default: True
    :return: read_counts_fil [str/file] preprocessed profile

    '''
    0 if os.path.exists(outdir) else os.mkdir(outdir)
    read_counts_fil = os.path.join(outdir, prefix + '_profile.xls')
    
    safile_fg, safile_bg = convert_saf(bedfil, fg=True), convert_saf(bedfil, fg=False)
    cntfiles, bams_join  = [], ' '.join(bamfils)
    
    for saf in [safile_fg, safile_bg]:
        if not bg and saf == safile_bg: break

        label = 'foreground' if saf == safile_fg else 'background'
        LOGS.info('Computing {} counts of peaks for all samples......'.format(label))
        cntfile = create_tmp_files('_featureCounts.cnt', tmpdir=kargs.tmpdir)[0].name
        
        cmd = 'featureCounts -F SAF -a {} -Q {} -T {} -o {} {} -O '.format(
                saf, kargs.mapq, kargs.thread, cntfile, bams_join)
        _sp = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = _sp.communicate()
        cntfiles.append(cntfile)
    
    LOGS.info('Writing pure samples profile into {}'.format(read_counts_fil))
    readcounts_matrix(cntfiles, kargs.infos['CELL'].values.tolist(), read_counts_fil)
    return read_counts_fil

def readcounts_matrix(count_files, cells, outfile, bak_offset=500000):
    '''
    convert read counts of peaks to matrix for total samples
    :param count_files: [str/files] read counts file of all samples
    :param cells: [list] cell or mixed sample names 
    :param outfile: [str/file path] output file
    :param bak_offset: [int] offset distance around summit to compute background read counts, default: 500000
    :return: 0
    
    '''
    count_files = count_files if isinstance(count_files, list) else [ count_files ]
    fps = [ open(fil, 'rb') for fil in count_files ]
    [ fp.readline() for fp in fps for idx in range(2) ]
    with open(outfile, 'w') as fp_w:
        fp_w.write('chrom\tstart\tend\t{}\n'.format('\t'.join(cells)))
        if len(count_files) > 1:
            while 1:
                for_cnts, bak_cnts = [fp.readline().split('\t') for fp in fps]
                if not for_cnts[0]: break
                tmp_infos  = for_cnts[1 : 4]
                cnt1, cnt2 = map(int, for_cnts[6 : ]), map(int, bak_cnts[6 : ])
                sig_ratio  = np.array(cnt1) / np.array(cnt2).astype(float) * bak_offset * 2 / 500.0
                sig_ratio  = np.around(sig_ratio, 3).tolist()
                tmp_infos += map(str, sig_ratio)
                fp_w.write('\t'.join(tmp_infos) + '\n')
        else:
            for idx, line in enumerate(fps[0]):
                line_infos = line.decode().strip().split('\t')
                tmp_infos = line_infos[1 : 4] + line_infos[6 : ]
                fp_w.write('\t'.join(tmp_infos) + '\n')
    return 0

def index_bamfile(bamfiles):
    '''
    check if the index file exists, otherwise build index
    :param bamfiles: [str/list] bam files
    ed_counts_filretrurn: 0
    
    '''
    bamfiles = bamfiles if isinstance(bamfiles, list) else [ bamfiles ]
    for bam in bamfiles:
        bam_index = bam + '.bai'
        if os.path.exists(bam_index): continue
        LOGS.warn('{} index cannot be detected, create index...'.format(bamfile))
        pysam.index(bam, bam_index) # create bam index
    return 0

def get_readcounts(bamfiles, paired=True):
    '''
    get total read counts of the bam
    :param bamfiles: [list/files] bam file pathes
    :param paired: [bool] only count paire-end reads, default: True
    :return: counts [int] total counts
    
    '''
    index_bamfile(bamfiles)
    counts = [ pysam.AlignmentFile(bam, 'rb').count() for bam in bamfiles ] # counts total reads of the bam files
    return counts
