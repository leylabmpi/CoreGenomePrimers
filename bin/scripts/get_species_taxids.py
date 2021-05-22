#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
from distutils.spawn import find_executable
from subprocess import Popen, PIPE

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'get species-level taxids'
epi = """DESCRIPTION:
For each provided taxid:
 * species-level taxids obtained via get_species_taxids.sh
All species-level taxids written to STDOUT (one taxid per line).
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('taxid_file', metavar='taxid_file', type=str,
                    help='One taxid per line')
parser.add_argument('--version', action='version', version='0.0.1')

def load_taxids(infile):
    taxids = set()
    with open(infile) as inF:
        for line in inF:
            line = line.strip()
            if line == '':
                continue
            taxids.add(line)
    return list(taxids)

def call_get_species_taxids(taxid):
    logging.info('Getting species-level taxids for {}'.format(taxid))
    p = Popen(['get_species_taxids.sh', '-t', str(taxid)], stdout=PIPE)
    output, err = p.communicate()
    rc = p.returncode

    taxids = list()
    for line in output.decode().split('\n'):
        line = line.strip()
        if line == '':
            continue
        taxids.append(line)
    logging.info('  No. of spec-level taxids: {}'.format(len(taxids)))
    return taxids
    
def get_spec_taxids(taxids):
    spec_taxids = set()
    for taxid in taxids:
           for x in call_get_species_taxids(taxid):
               spec_taxids.add(x)
    return spec_taxids

def main(args):
    # checking that executables exist
    exe = ['blastdbcmd', 'get_species_taxids.sh']
    for x in exe:
        if find_executable(x) is None:
            raise OSError('Cannot find executable: {}'.format(x))
    # load taxids
    taxids = load_taxids(args.taxid_file)
    # get species-level taxids
    taxids = get_spec_taxids(taxids)
    # writing output
    for taxid in taxids:
        print(taxid)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
