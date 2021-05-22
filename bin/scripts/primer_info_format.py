#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
from collections import defaultdict

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Template python script'
epi = """DESCRIPTION:
* placeholder1
* placehodler2
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('primer_info_file', metavar='primer_info_file', type=str,
                    help='Primer info table file')
parser.add_argument('--version', action='version', version='0.0.1')


def _open(infile, mode='rb'):
    """
    Openning of input, regardless of compression
    """
    if infile.endswith('.bz2'):
        return bz2.open(infile, mode)
    elif infile.endswith('.gz'):
        return gzip.open(infile, mode)
    else:
        return open(infile)

def _decode(x):
    """
    Decoding input, if needed
    """
    try:
        x = x.decode('utf-8')
    except AttributeError:
        pass
    return x

def parse_table(infile):
    header = {}
    primers = defaultdict(dict)
    with open(infile) as inF:
        for i,line in enumerate(inF):
            line = line.rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                continue
            primer_set = line[header['primer_set']]
            primer_type = line[header['primer_type']]
            primer_seq = line[header['sequence']]
            primers[primer_set][primer_type] = primer_seq
    return(primers)

def write_primer_table(primers):
    types = ['PRIMER_LEFT', 'PRIMER_RIGHT']
    for primer_set in primers.keys():
        seqs = []
        for T in types:
            try:
                seqs.append(primers[primer_set][T])
            except KeyError:
                continue
        if len(seqs) == 2:
            print('\t'.join([str(primer_set)] + seqs))

def main(args):
    primers = parse_table(args.primer_info_file)
    write_primer_table(primers)
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
