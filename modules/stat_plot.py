#!/usr/bin/env python
#title      : stat_plot.py
#decription : Visualize the results of deconATAC
#author     : Huamei Li
#date       : 23/07/2018
#type       : module
#version    : 3.8

#-----------------------------------------------------
# load python modules

import sys
import matplotlib as mpl
if sys.platform == 'win32':
	pass
else:
	mpl.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import cycle
from rpy2.robjects import r, pandas2ri

#-----------------------------------------------------
# load own modules

from modules.utils import *

#-----------------------------------------------------
# global setting

plt.rcParams["font.family"] = "Times New Roman"
sns.set(color_codes=True)

#-----------------------------------------------------

def cluster_heatmap(df, outfile):
    '''
    plot cluster heatmap
    :param df: [pd.DataFrame] dataframe of input data
    :param outfile: [str] output graph file name 
    :return: 0
    
    '''
    outfile = outfile if outfile.endswith('.png') else outfile + '.png'
    r.assign('outfile', outfile)
    r.assign('df', df)

    r('''
        options(warn = -1)
        list.of.packages <- c("colorRamps")
        new.pkgs <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
        if(length(new.pkgs )) install.packages(new.pkgs)
        suppressMessages(library(colorRamps, quietly = T))
        df = as.matrix(df)
        png(outfile)
        heatmap(df, col = blue2red(50), margins=c(15, 5), labRow=NA)
        dev.off()
    ''')
    
    return 0

def stack_bars(df, outfile, tool='deconPeaker'):
    '''
    stack bar plot for mixture cell popuplation
    :param df: [pd.DataFrame] rows: cell proportation; columns: mixed samples
    :param outfile: [str] output file name
    :param tool: [str] deconvolution tool, default: deconATAC 
    :return: 0
    
    '''
    nrows, ncols = df.shape
    colors = plt.cm.Set1(np.linspace(0, 1, ncols))
    df.plot(kind='bar', stacked=True, rot=45, figsize=(10, 6), color=colors)
    
    plt.ylim([0, 1])
    plt.xlabel('Sequential substraction of cell types', fontsize=10)
    plt.ylabel(tool + ' prediction', fontsize=10)
    plt.yticks(np.arange(0, 1.2, step=0.2), ('0', '0.2', '0.4', '0.6', '0.8', '1.0'))
    plt.tight_layout() 
    outfile = outfile if outfile.endswith('.png') else outfile + '.png'
    plt.savefig(outfile, format='png', dpi=300)
    return 0

def bars(df, outfile, platform='ATAC-Seq'):
    '''
    bar plot for the pandas dataframe
    :param df: [pd.DataFrame] data which need to be plot
    :param outfile: [str/file] output file name
    :return: 0
    
    '''
    df.plot.bar(x='lab', y='val', rot=45, legend=False, fontsize=10)
    plt.xlabel('')
    plt.tight_layout()
    if platform == 'ATAC-Seq':
        plt.title('Number of cell type specific peaks', fontsize=12)
    elif platform == 'RNA-Seq':
        plt.title('Number of cell type specific genes', fontsize=12)
    else:
        plt.title('Number of cell type specific probes', fontsize=12)
    outfile = outfile if outfile.endswith('.png') else outfile + '.png'
    plt.tight_layout()
    plt.savefig(outfile, format='png', dpi=300)
