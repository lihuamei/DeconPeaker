#!/usr/bin/env python
#title      : stat_plot.py
#decription : Visualize the results of deconATAC
#author     : Huamei Li
#date       : 23/07/2018
#type       : module
#version    : 2.7

#-----------------------------------------------------
# load python modules

import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import cycle
from rpy2.robjects import r, pandas2ri

#-----------------------------------------------------
# load own modules

from utils import *

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
        library(colorRamps)
        df = as.matrix(df)
        #row_order = hclust(dist(df))$order
        #col_order = hclust(dist(t(df)))$order
        #df_hclust = df[row_order, col_order]
        png(outfile)
        heatmap(df, col = blue2red(50), margins=c(15, 5), labRow=NA)
        dev.off()
    ''')
    
    return 0

def scatter_true_est(df, rsquare, labels=None, outfile=None):
    '''
    scatter plot for true and estimate population
    :param df: [pd.DataFrame] dataframe of input data, first column: True population, second colum: estimate population 
    :param rsquare: [float] rsqaure index of prediction
    :param labels: [list] labels of x and y axes, default: None
    :param outfile: [str/file path] output graph file name
    
    '''
    if not labels: axis_labels=['Synthetic ground truth', 'deconATAC prediction']
    markers = ['.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '*', 'h', '+', 'H']
    colors  = ['tan', 'darkgreen', 'deepskyblue', 'navy', 'plum', \
            'olive', 'crimson', 'green', 'lightpink', 'orchid', 'teal', 'yellow']
    marker_gens, color_gen = cycle(markers), cycle(colors)
    (nrows, ncols), labels, columns = df.shape, df.index.values, df.columns
    colors, markers = [], []
    
    [(colors.append(next(color_gen)), markers.append(next(marker_gens))) for idx in xrange(nrows) ]
    
    fig = plt.figure(figsize=(10, 6))
    for index, row in enumerate(df.values):
        plt.scatter(row[0], row[1], c=colors[index], s=10, marker=markers[index], label=labels[index])

    plt.xlim([0, 1]); plt.ylim([0, 1])
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=nrows/20 + 1)
    plt.plot(range(2), '--', alpha=0.5, linewidth=1)
    
    plt.xticks(np.arange(0, 1.2, step=0.2), ('0', '0.2', '0.4', '0.6', '0.8', '1.0'))
    plt.yticks(np.arange(0, 1.2, step=0.2), ('', '0.2', '0.4', '0.6', '0.8', '1.0'))
    if rsquare >= 0:
        plt.text(0.05, 0.85, '$R^2 = {} $'.format(rsquare), color='blue')

    plt.xlabel(axis_labels[0], fontsize=12)
    plt.ylabel(axis_labels[1], fontsize=12)
    plt.tight_layout() 

    outfile = outfile if outfile.endswith('.png') else outfile + '.png'
    plt.savefig(outfile, format='png', dpi=300)
    return 0

def stack_bars(df, outfile, tool='deconATAC'):
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
    leg = plt.legend(
            loc = 'center left'      , 
            bbox_to_anchor = (1, 0.5), 
            ncol = nrows/20 + 1      ,
            labelspacing = 0.3       ,
            borderpad = 0.5          , 
            handletextpad = 0.3      ,
            handlelength = 0.8       ,
            columnspacing = 0.8      ,
        )
    
    plt.xlabel('Sequential substraction of cell types', fontsize=10)
    plt.ylabel(tool + ' prediction', fontsize=10)
    plt.yticks(np.arange(0, 1.2, step=0.2), ('0', '0.2', '0.4', '0.6', '0.8', '1.0'))
    plt.gca().axes.xaxis.set_ticklabels([])
    plt.tight_layout() 
    
    outfile = outfile if outfile.endswith('.png') else outfile + '.png'
    plt.savefig(outfile, format='png', dpi=300)
    return 0

def bars(df, outfile):
    '''
    bar plot for the pandas dataframe
    :param df: [pd.DataFrame] data which need to be plot
    :param outfile: [str/file] output file name
    :return: 0
    
    '''
    df.plot.bar(x='lab', y='val', rot=45, legend=False, fontsize=10)
    plt.xlabel('')
    plt.tight_layout()
    plt.title('Number of cell type specific peaks', fontsize=12) 
    outfile = outfile if outfile.endswith('.png') else outfile + '.png'
    plt.savefig(outfile, format='png', dpi=300)

#if __name__ == '__main__':
    #iris = sns.load_dataset("iris")
    #species = iris.pop("species")
    #iris = iris.iloc[0:50]
    #iris = iris/np.sum(iris)
    #outfile = 'test.png'
    #data = pd.DataFrame({'val': [982, 745, 2093, 3962, 7984, 10042, 13151, 8410, 88, 8635, 6426, 885, 1940], 'lab': ['CD8T', 'GMP', 'CMP', 'CD4T', 'Ery', 'CLP', 'MONO', 'B', 'MPP', 'MEP', 'NK', 'HSC', 'LMPP']})
    #bars(data, outfile)
