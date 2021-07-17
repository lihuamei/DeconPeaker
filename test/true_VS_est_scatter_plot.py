#!/usr/bin/env python
#title       : true_VS_est_scatter_plot.py
#description : plot scatter based on true and estimate proportions.
#author      : Huamei Li
#date        : 19/06/2018
#type        : module
#version     : 2.7

#-----------------------------------------------------
# load python module

import os
import sys, pdb
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import cycle
import seaborn as sns; sns.set()
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
#-----------------------------------------------------
# global setting

#from matplotlib.font_manager import _rebuild; _rebuild()
#plt.rcParams["font.family"] = "Times New Roman"

#-----------------------------------------------------

if len(sys.argv[1 : ]) != 3:
    sys.exit('[ERROR] USAGE: python true_VS_est_scatter_plot.py TRUE_proportion_file Estimate_proportion_file Algorithm_Name')
else:
    print('[INFO] Plotting...')

true_prop = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
esti_prop = pd.read_csv(sys.argv[2], sep='\t', header=0, index_col=0)
prefix_na, colnames = sys.argv[-1], true_prop.columns & esti_prop.columns

esti_prop, true_prop = esti_prop[colnames], true_prop[colnames]
#esti_prop.index = [line.split('_', 1)[-1] for line in esti_prop.index]
rownames = true_prop.index & esti_prop.index
esti_prop, true_prop = esti_prop.loc[rownames], true_prop.loc[rownames]

corr, pvalue = pearsonr(true_prop.values.flatten(), esti_prop.values.flatten())
colors = np.random.rand(3, true_prop.shape[1])
for index, name in enumerate(colnames):
    x, y = true_prop[name], esti_prop[name]
    plt.scatter(x = x, y = y, s=10, c=colors[:, index], label=name)

plt.xlim([0, 1]); plt.ylim([0, 1])
plt.xticks(np.arange(0, 1.2, step=0.2), ('0', '0.2', '0.4', '0.6', '0.8', '1.0'))
plt.yticks(np.arange(0, 1.2, step=0.2), ('', '0.2', '0.4', '0.6', '0.8', '1.0'))
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.plot(range(2), '--', alpha=0.5, linewidth=1)

axis_labels=['Truth', 'Prediction']
plt.xlabel(axis_labels[0], fontsize=12)
plt.ylabel(axis_labels[1], fontsize=12)

plt.text(.015, 0.90, r'p.value={}'.format(pvalue), fontsize=8, color='b', fontweight='bold')
plt.text(.015, 0.95, r'correration={}'.format(corr), fontsize=8, color='b', fontweight='bold')
plt.text(.015, 0.85, r'n={}'.format(esti_prop.shape[0]), fontsize=8, color='b', fontweight='bold')
plt.tight_layout()
plt.savefig(prefix_na + '_true_estimate_scatter.png', dpi=300, format='png')
