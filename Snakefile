# import 
from __future__ import print_function
import os
import sys
import socket
import getpass
import random
import numpy as np
import pandas as pd
from distutils.spawn import find_executable

# --setup--
## load
configfile: 'config.yaml'
## pipeline utils
snake_dir = config['pipeline']['snakemake_folder']
include: snake_dir + 'bin/ll_pipeline_utils/Snakefile'
config_default(config, 'pipeline', 'name')
## output dir
config['output_dir'] = config['output_dir'].rstrip('/') + '/'
## samples table
if not os.path.isfile(config['samples_file']):
    raise IOError('Cannot find file: {}'.format(config['samples_file']))
config['samples'] = pd.read_csv(config['samples_file'], sep='\t', comment='#')
for x in ['Taxon', 'Taxid', 'Fasta', 'Domain']:
    if x not in config['samples'].columns:
        msg ='Cannot find "{}" column in the samples table'
        raise ValueError(msg.format(x))
### Fasta files
for F in config['samples'].Fasta.tolist():
    if not os.path.isfile(F):
        msg = 'Cannot find file: {}'
        raise IOError(msg.format(F))
### samples
config['samples'].Taxon = config['samples'].Taxon.astype(str)
config['samples'] = config['samples'].set_index(config['samples'].Taxon)
config['samples_unique'] = config['samples'].Taxon.unique().tolist()
### taxids
config['taxids'] = config['samples'].Taxid.unique().tolist()
### filters
if 'Filters' in config['samples'].columns:
    config['filters'] = '--filter-list ' + ','.join(config['samples'].Filters.unique().tolist())
else:
    config['filters'] = ''

## temp directory
config['pipeline']['username'] = getpass.getuser()
config['pipeline']['email'] = config['email']
config['tmp_dir'] = os.path.join(config['tmp_dir'], config['pipeline']['username'])
config['tmp_dir'] = os.path.join(config['tmp_dir'], 'CoreGenomePrimers_' + str(os.stat('.').st_ino) + '/')
print('\33[33mUsing temporary directory: {} \x1b[0m'.format(config['tmp_dir']))
print('\33[33mUsing output directory: {} \x1b[0m'.format(config['output_dir']))
if not os.path.isdir(config['tmp_dir']):
    os.makedirs(config['tmp_dir'])
    if not os.path.isdir(config['tmp_dir']):
        raise IOError('Cannot find tmp dir: {}'.format(config['tmp_dir']))

## including modular snakefiles
include: snake_dir + 'bin/dirs'
include: snake_dir + 'bin/Snakefile'
include: snake_dir + 'bin/utils/Snakefile'
include: snake_dir + 'bin/taxids/Snakefile'
include: snake_dir + 'bin/cgp/genes/Snakefile'
include: snake_dir + 'bin/cgp/clusters/Snakefile'
include: snake_dir + 'bin/cgp/nontarget/cds/Snakefile'
include: snake_dir + 'bin/cgp/nontarget/rrna/Snakefile'
include: snake_dir + 'bin/cgp/primers/Snakefile'
include: snake_dir + 'bin/cgp/amplicon/Snakefile'

## local rules
localrules: all

# rules
rule all:
    input:
        which_rules

