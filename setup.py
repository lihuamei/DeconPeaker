#!/usr/bin/env python

import os
from setuptools import setup, find_packages

def before_install():
    cur_dir = os.getcwd()
    os.system('export PYTHONPATH=%s:$PYTHONPATH'.format(cur_dir))

before_install()

setup(
    name = 'DeconPeaker',
    version = '1.1.1',
    packages = find_packages(),
    install_requires=[
        'bx==0.3.0',
        'matplotlib==3.7.5',
        'bx-python',
        'numexpr',
        'numpy',
        'pandas',
        'pysam',
        'rpy2',
        'scipy',
        'seaborn',
        'setuptools'
    ],
    entry_points = {
        'console_scripts': [
            'deconPeaker=DeconPeaker.deconPeaker:run',
        ],
    },
    package_data = {
        '': ['data/*', 'extra/*', 'modules/*', '*.r']
    },
    author = 'Huamei Li',
    author_email = 'li_hua_mei@163.com',
    description = 'A deconvolution method to estimate cell type proportions in chromatin accessibility data',
    url = 'https://github.com/lihuamei/DeconPeaker.git',
    include_package_data = True
)
