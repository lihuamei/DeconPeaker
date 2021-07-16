#!/usr/bin/env python
#title       : random_proportion.py 
#description : Generate random proportion profile for simulating mixture samples. 
#author      : Huamei Li
#date        : 19/09/2018
#type        : module
#version     : 3.8

#-----------------------------------------------------
# load own modules

from modules.utils import *

#-----------------------------------------------------

def random_proportions(sndf, mixture, rep_counts, prefix, outdir):
    '''
    Generate random proportions for mixture samples
    :param sndf: [pd.DataFrame] pure cell informative, include file path, cell name
    :param mixture: [list] number of cells need to be mixed
    :param rep_counts: [int] number of mixed samples randomly generated from fixed cell types
    :param prefix: [str] prefix name of output files
    :param outdir: [str] output directory
    :return: random_results 
    
    '''
    np.random.seed(100)
    cellnum, cellnames = sndf.shape[0], sndf.cellname
    random_nums  = mixture if mixture else range(1, cellnum + 1)
    random_cells = [ np.random.choice(cellnum, num, replace=False) \
            for num in random_nums for idx in xrange(rep_counts) ] 
    random_props   = [np.random.dirichlet([1] * len(lst), 1) for lst in random_cells]
    random_results = pd.DataFrame(np.zeros((len(random_cells),cellnum)), columns=sndf.cellname)
    
    for idx, cells in enumerate(random_cells):
        props = random_props[idx][0]
        for cell, prop in zip(cells, props):
            random_results.loc[idx, cellnames[cell]] = prop
    
    sample_labels = [ 'Sample_{}'.format(idx + 1) for idx in xrange(len(random_cells)) ]
    random_results.index = sample_labels
    filename = os.path.join(outdir, prefix + '_Dirichlet_proportions.xls')
    random_results.to_csv(filename, sep='\t', header=True, index=True)
    return random_results

