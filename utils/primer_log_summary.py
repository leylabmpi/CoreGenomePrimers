#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
from pathlib import Path

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Summarize info from primer design/filter log files'
epi = """DESCRIPTION:
Quickly create a tab-delim table of log file info for
the primer design and filtering steps in the pipeline.

The input should be the base directory to the pipeline output
logs ($OUTDIR/logs/).

The output table is written to STDOUT.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('base_dir', metavar='base_dir', type=str,
                    help='base directory for log files')
parser.add_argument('--version', action='version', version='0.0.1')


def parse_log(path, stage):    
    gene = path.parts[-2]
    cluster = os.path.splitext(path.name)[0]
    with open(str(path)) as inF:
        for line in inF:
            line = line.rstrip().split(' - ', 2)
            if (line[1].lstrip().startswith('File written:') or
                line[1].startswith('Reading file:') or
                line[1].lstrip().startswith('No. of records')):
                continue
            line = [stage, gene, cluster, line[1]]
            print('\t'.join(line))

def main(args):
    paths = list(Path(args.base_dir).rglob('*.log'))
    print('\t'.join(['stage', 'gene_type', 'cluster', 'log']))
    # design    
    for p in paths:        
        if p.parts[-3] == 'design':     
            parse_log(p, 'design')
    # crosshyb
    for p in paths:
        if p.parts[-3] == 'blastn_crosshyb_filter':     
            parse_log(p, 'crosshyb')
    # non-target
    for p in paths:
        if p.parts[-3] == 'blastn_other_taxa_filter':     
            parse_log(p, 'non-target')
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
