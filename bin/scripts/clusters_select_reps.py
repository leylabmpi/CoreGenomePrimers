#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Selecting representative sequences'
epi = """DESCRIPTION:
* placeholder1
* placehodler2
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('reps_tsv', metavar='reps_tsv', type=str,
                    help='mmseqs cluster reps as tsv')
parser.add_argument('genes_fasta', metavar='genes_fasta', type=str,
                    help='genes nucleotide fasta')
parser.add_argument('--version', action='version', version='0.0.1')



def main(args):
    # reps
    logging.info('Reading file: {}'.format(args.reps_tsv))
    reps = set()
    with open(args.reps_tsv) as inF:
        for line in inF:
            if line.rstrip() == '':
                continue
            reps.add(line.split('\t')[0])
    logging.info('  No. of rep. seqs.: {}'.format(len(reps)))
    # filtering fasta
    logging.info('Reading file: {}'.format(args.genes_fasta))
    cnt = 0
    with open(args.genes_fasta) as inF:
        to_keep = False
        for line in inF:
            line = line.rstrip()
            if line.startswith('>'):
                if line.lstrip('>') in reps:
                    to_keep = True
                    cnt += 1
                else:
                    to_keep = False
            if to_keep is True:
                print(line)
    if cnt == 0:
        raise ValueError('No rep. sequences found!')
    logging.info('  No. of seqs. written: {}'.format(cnt))
                
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
