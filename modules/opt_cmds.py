#!/usr/bin/env python
#title       : opt_cmd.py
#description : Required parameters for program operation
#author      : Huamei Li
#date        : 29/05/2018
#type        : module
#version     : 3.8

#--------------------------------------------------------
# load python modules

import sys
import argparse

#--------------------------------------------------------
# define options for deconPeaker

def opts():
    parser = argparse.ArgumentParser(description='deconPeaker - a deconvolution model to identify cell types based on \
            chromatin accessibility in ATAC-Seq data of mixture samples.')
    subpar = parser.add_subparsers() # add sub parsers

    preprocess = subpar.add_parser('preprocess' , help='Create chromatin accessibility profile for pure samples, \
            this step is the basis for subsequent specific cell type identification and mixture deconvolution. Note: This step only support Linux system, \
            and if you have a large sample size for pure cells, please ensure enough sufficient memory and hard storage space for program to run normally.')
    
    preprocess.add_argument(
            '--pure', 
            '-p', 
            help = 'YAML format file, details about a set of pure samples, including the pathes of BAM and peak files (MACS2). Please refer \
                    to the test configure YAML of pure samples: test/test_pure_samples.yaml',
            type = str,
            metavar = 'PURE'
        )

    preprocess.add_argument(
            '--hg-version',
            '-g',
            help = 'The genomic version specified by the sequence alignment. DEFAULT: hg19',
            choices = ['hg19', 'hg38'],
            default = 'hg19'
        )

    preprocess.add_argument(
            '--qvalue',
            help = 'Filtering weak peaks by qvalue, meaning that peaks above this value need to be filtered out. DEFAULT: 0.01',
            type = float,
            metavar = 'QVALUE',
            default = 0.01
        )
    
    preprocess.add_argument(
            '--thread',
            '-t',
            help = 'Number of thread. DEFAULT: 4',
            type = int,
            default = 4
        ) 
           
    preprocess.add_argument(
            '--mapq',
            '-q',
            help = 'The minimum mapping quality score a read must satisfy in order to be counted. \
                    For paired-end reads, at least one end should satisfy this criteria. DEFAULT: 30',
            type = int,
            metavar = 'MAPQ',
            default = 30
        )
    
    preprocess.add_argument(
            '--offset',
            '-s',
            help = 'Offset distance from TSS. Note: peaks falling in the area [TSS - offset, TSS + offset] will be removed. DEFAULT: 1500 bp',
            type = int,
            metavar = 'OFFSET',
            default = 1500
       )

    preprocess.add_argument(
            '--prefix',
            help = 'Prefix name of preprocessed results. DEFAULT: None',
            type = str,
            metavar = 'PREFIX',
       )
    
    preprocess.add_argument(
            '--outdir',
            '-o',
            help = 'If specified all output files will be written to that directory. DEFAULT: the current working directory',
            type = str,
            metavar = 'OUTDIR',
            default = './'
        )

    findctsps = subpar.add_parser('findctsps', help='Find cell type specific peaks/genes accross pure samples, different pure cell samples \
            require replicates as input.')

    findctsps.add_argument(
            '--lib-strategy',
            '-l',
            help = 'A string indicating the type of the profile measurements. Used to choose profile loading method and fine-tuning of \
                    analytical methods. DEFAULT: ATAC-Seq',
            choices = ['ATAC-Seq', 'RNA-Seq', 'Microarray'],
            default = 'ATAC-Seq'
        )
    
    findctsps.add_argument(
            '--profile',
            '-f',
            help = 'Profile of peaks in pure cell samples,  row = peaks, column = pure cell names, \
                    each cell is accessibility of peaks. NOTE: ATAC-Seq data input format, first three columns are chromosome, \
                    start and end position respectively; RNA-Seq and Microarray: first column must be the genes/probes',
            type = str,
            metavar = 'PROFILE'
        )

    findctsps.add_argument(
            '--phenotype',
            '-c',
            help = 'Phenotype classes file, which columns correspong in exact order to the reference samples in the \
                    reference samples file and rows correspond to the cell type classes taht will be used to definethe \
                    classes in the cell specific peaks file. value 1 = indicate membership of the reference sample, \
                    value 2 = indicates the class that the sample will be compared against.',
            type = str,
            metavar = 'PHENOTYPE'
        )
    
    findctsps.add_argument(
            '--norm',
            help = 'A series method for normalizing pure samples read counts, including quantile, DESeq, UQN (upper quantile), PPM and TMM. \
                    DEFAULT: QN',
            choices = ['QN', 'DESeq', 'UQN', 'PPM', 'TMM'],
            default = 'QN'
        )

    findctsps.add_argument(
            '--min-group-size',
            '-m'              ,
            help = 'Minimum number of cell type specific peaks to consider from each phenotype for signature matrix. DEFAULT: 100',
            metavar = 'MIN-GROUP-SIZE',
            type = int, 
            default = 100
        )

    findctsps.add_argument(
            '--max-group-size',
            '-x'              ,
            help = 'Maximum number of cell type specific peaks to consider from each phenotype for signature matrix. DEFAULT: 200',
            metavar = 'MAX-GROUP-SIZE',
            type = int,
            default = 200
        )

    findctsps.add_argument(
            '--score',
            '-s'      ,
            help = 'Pi-score value for high varibale peaks. DEFAULT: 1',
            metavar = 'SCORE',
            type = float,
            default = 1
        )

    findctsps.add_argument(
            '--ratio',
            help = 'Peak expression concentration ratio of top two cell types with the highest expression of \
                    each peak, the ratio of the second highest to the highest is less than 0.33, and will be \
                    filtered out. DEFAULT: 0.33',
            metavar = 'RATIO',
            type = float,
            default = 0.33
        )

    findctsps.add_argument(
            '--thread',
            '-t',
            help = 'Number of thread. DEFAULT: 4',
            type = int,
            default = 1
        )

    findctsps.add_argument(
            '--outdir',
            '-o',
            help = 'If specified all output files will be written to that directory. DEFAULT: the current working directory',
            type = str,
            metavar = 'OUTDIR',
            default = './'
        )
    
    deconv = subpar.add_parser('deconvolution', help='Based on pure cell profile information, robust regression deconvolution \
            strategy was used to estimate the proportion of possible cell types in the mixed samples.')
    deconv.add_argument(
            '--lib-strategy',
            '-l',
            help = 'A string indicating the type of the profile measurements. Used to choose profile loading method and fine-tuning of \
                    analytical methods. DEFAULT: ATAC-Seq',
            choices = ['ATAC-Seq', 'RNA-Seq', 'Microarray'],
            default = 'ATAC-Seq'
        )
    
    deconv.add_argument(
            '--mixture',
            '-m',
            help = 'Support BAM and TABLE format file. [TABLE] Mixture file (GEP matrix: row 1 = sample labels; \
                    column 1 = gene symbols; no missing values), [BAM] Generated by alignment tools, multiple bams \
                    must be saved in YANL format file. Note: reference YAML -> test/mixture_samples.yaml',
            type = str,
            metavar = 'MIXTURE',
        )

    deconv.add_argument(
            '--pure',
            '-p',
            help = 'Signature matrix file (Cell type PEP matrix: row 1 = sample labels; column 1 = peaks / gene symbols; no missing values)',
            type = str,
            metavar = 'PURE'
        )
   
    deconv.add_argument(
            '--format',
            '-f',
            help = 'Format of mixture files, support TABLE and BAM formats. DEFAULT: TABLE',
            choices = ['TABLE', 'BAM'],
            default = 'TABLE'
        )

    deconv.add_argument(
            '--method',
            help = 'Deconvolution method, including SIMPLS, LR (Linear regression) and RSIMPLS (Robust SIMPLS). DEFAULT: SIMPLS',
            choices = ['SIMPLS', 'RSIMPLS', 'LR'],
            default = 'SIMPLS'
        )

    deconv.add_argument(
            '--pvalue',
            help = 'Estimate P-value of deconvolution or not. DEFAULT: FALSE',
            choices = ['TRUE', 'FALSE'],
            default = 'FALSE'
        )

    deconv.add_argument(
            '--norm',
            help = 'A series method for normalizing mixture samples, including quantile, DESeq, upper quantile, PPM and TMM. \
                    DEFAULT: None',
            choices = ['QN', 'DESeq', 'UQN', 'PPM', 'TMM'],
            default = None
        )

    deconv.add_argument(
            '--thread',
            '-t',
            help = 'Number of thread for counting reads of BAMs. DEFAULT: 1',
            type = int,
            default = 1
        )

    deconv.add_argument(
            '--outdir',
            '-o',
            help = 'If specified all output files will be written to that directory. Default: the current working directory',
            type = str,
            metavar = 'OUTDIR',
            default = './'
        )

    simulate = subpar.add_parser('simulation', help='Simulate mixed samples with different proportions of cells. \
            [Note] The method is proportional random sampling of reads from different cell types of BAM/BED files.')
    
    simulate.add_argument(
            '--pure',
            '-p',
            help = 'YAML format file, details about a set of pure samples, including the pathes of BAMS. \
                    Note: reference YAML -> test/synthetic_mixture_samples.yaml',
            type = str,
            metavar = 'PURE'
        )

    simulate.add_argument(
            '--mixture',
            '-m',
            help = 'Number of mixed cell types. Note: Generating a number of samples of different cell types requires separating \
                    the different numbers by commas (such as 3,4,5). If the number of cell types in the sample is not specified, \
                    the program will randomly generate all possiable mixed samples. DEFAULT: None',
            type = str,
            metavar = 'MIXTURE',
            default = None
        )

    simulate.add_argument(
            '--rep-counts',
            '-r',
            help = 'Repeat generating the number of random samples containing a specific combination of cell types. DEFAULT: 5',
            metavar = 'REPLICATES',
            type = int,
            default = 5
        )

    simulate.add_argument(
            '--readcounts',
            '-c',
            help = 'Total read counts of each simulated mixture sample. DEFAULT: 2000000',
            metavar = 'READCOUNTS',
            type = int,
            default = 2000000
        )

    simulate.add_argument(
            '--format',
            '-f',
            help = 'Format of tag file of pure cells, support "BAM" and "BED". NOTE: make sure samtools and bedtools are installed',
            choices = ['BAM', 'BED'],
        )

    simulate.add_argument(
            '--genome',
            '-g'      ,
            help = 'Reference genomic file location, only corresponding with --format=BED.',
            metavar = 'GENOME',
            type = str
        )

    simulate.add_argument(
            '--thread',
            '-t',
            help = 'Number of thread. DEFAULT: 4',
            type = int,
            default = 4
        )

    simulate.add_argument(
            '--prefix',
            help = 'Prefix name of simulated result files, default: deconPeaker_Sim',
            metavar = 'PREFIX',
            type = str,
            default = 'deconPeaker_Sim'
        )

    simulate.add_argument(
            '--outdir',
            '-o',
            help = 'If specified all output files will be written to that directory. DEFAULT: the current working directory',
            type = str,
            metavar = 'OUTDIR',
            default = './'
        )


    if len(sys.argv) <= 1: sys.exit(parser.print_help())
    args = parser.parse_args()
    if len(sys.argv) == 2:
        if sys.argv[-1] == 'preprocess': sys.exit(preprocess.print_help())
        if sys.argv[-1] == 'findctsps': sys.exit(findctsps.print_help())
        if sys.argv[-1] == 'deconvolution': sys.exit(deconv.print_help())
        if sys.argv[-1] == 'simulation': sys.exit(simulate.print_help())
    
    args.sub_parser, args.preprocess, args.findctsps, args.deconv, args.simulate = \
            sys.argv[1], preprocess, findctsps, deconv, simulate
    return args

if __name__ == '__main__':
    opts()
