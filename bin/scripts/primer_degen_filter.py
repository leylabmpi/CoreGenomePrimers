#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import shutil
import logging
import itertools
import collections
import functools
from subprocess import Popen, PIPE
from distutils.spawn import find_executable
# 3rd party
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

# arg parse
desc = 'Filter out undesireable primers'
epi = """DESCRIPTION:
Filtering degenerate primers based on parameters (eg., Tm) 
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('primer_degen_fasta', metavar='primer_degen_fasta', type=str,
                    help='Degenerate primer fasta file')
parser.add_argument('primer_expand_fasta', metavar='primer_expand_fasta', type=str,
                    help='Expanded primer sequence fasta file')
parser.add_argument('primer_info', metavar='primer_fino', type=str,
                    help='Primer info table')
parser.add_argument('--prefix', type=str, default='primers_filtered',
                    help='Output file prefix')
parser.add_argument('--n-primers-max', type=int, default=0,
                    help='Max number of degenerate primers kept (0 = all kept)')
parser.add_argument('--degen-max', type=int, default=64,
                    help='Max no. of degeneracies in the primer')
parser.add_argument('--Tm-min', type=float, default=55,
                    help='Min melting temp of either primer in the primer pair')
parser.add_argument('--Tm-max', type=float, default=65,
                    help='Max melting temp of either primer in the primer pair')
parser.add_argument('--Na', type=float, default=50,
                    help='Concentration of Na ions in the PCR [mM]')
parser.add_argument('--Tris', type=float, default=0,
                    help='Concentration of Tris in the PCR [mM]')
parser.add_argument('--Mg', type=float, default=0,
                    help='Concentration Mg ions in the PCR [mM]')
parser.add_argument('--dNTPs', type=float, default=0.2,
                    help='Concentration dNTPs in the PCR [mM]')
parser.add_argument('--version', action='version', version='0.0.1')



def calc_Tm(seq, Na=50, Tris=0, Mg=0, dNTPs=0):
    """
    Calculating melting temperature
    """
    return mt.Tm_NN(Seq(seq), Na=Na, Tris=Tris, Mg=Mg, dNTPs=dNTPs)

def filter_primer_info(infile, prefix, Tm_min=55, Tm_max=65, Na=50, Tris=0, Mg=0, dNTPs=0):
    """
    Filtering primers based on Tm
    """
    logging.info('Reading file: {}'.format(infile))
    outfile = prefix + '.tsv'
    bad_primers = set()
    status = {'total' : 0, 'low Tm' : 0, 'high Tm' : 0}
    with open(infile) as inF, open(outfile, 'w') as outF:
        header = {}
        for i,line in enumerate(inF):
            line = line.rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                line += ['Tm_fwd', 'Tm_rev']
            elif line[0] == '':
                continue
            else:
                status['total'] += 1
                primer_id = line[header['primer_id']]
                fwd_seq = line[header['F_Primer_IUPAC_Sequence']]
                rev_seq = line[header['R_Primer_IUPAC_Sequence']]
                # fwd primer
                Tm_fwd = calc_Tm(fwd_seq, Na, Tris, Mg, dNTPs)
                if Tm_fwd < Tm_min or Tm_fwd > Tm_max:
                    if Tm_fwd < Tm_min:
                        status['low Tm'] += 1
                    elif Tm_fwd > Tm_max:
                        status['high Tm'] += 1
                    bad_primers.add(primer_id)
                    continue
                # rev primer
                Tm_rev = calc_Tm(rev_seq, Na, Tris, Mg, dNTPs)
                if Tm_rev < Tm_min or Tm_rev > Tm_max:
                    if Tm_rev < Tm_min:
                        status['low Tm'] += 1
                    elif Tm_rev > Tm_max:
                        status['high Tm'] += 1
                    bad_primers.add(primer_id)
                    continue
                # summary
                line += [str(round(Tm_fwd,1)), str(round(Tm_rev,1))]
            outF.write('\t'.join(line) + '\n')
    # status
    msg = '  No. of primer pairs assessed: {}'
    logging.info(msg.format(status['total']))
    msg = '  No. of primer pairs with a non-desirable Tm: {}'
    logging.info(msg.format(len(bad_primers)))
    logging.info('    No. with a low Tm: {}'.format(status['low Tm']))
    logging.info('    No. with a high Tm: {}'.format(status['high Tm']))
    logging.info('File written: {}'.format(outfile))
    return bad_primers

def edit_primers_fasta(infile, prefix, filetype='degen', to_filter=set()):
    """
    Reading in primer fasta file
    """    
    logging.info('Reading file: {}'.format(infile))
    outfile = prefix + '_{}.fna'.format(filetype)
    if filetype == 'degen':
        regex = re.compile(r'([0-9]+)([fr])')
    else:
        regex = re.compile(r'([0-9]+)[fr].*')
    seq_id,fr,rep_id = None,None,None
    to_keep = False
    cnt = 0
    with open(infile) as inF, open(outfile, 'w') as outF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            if line.startswith('>'):
                seq_id = regex.search(line).groups()[0]
                if seq_id in to_filter:
                    to_keep = False
                else:
                    to_keep = True
                    cnt += 1
            if to_keep is True:
                outF.write(line + '\n')                             
    logging.info('  File written: {}'.format(outfile))
    logging.info('    No. of seqs written: {}'.format(cnt))

def main(args):
    outdir = os.path.split(args.prefix)[0]
    if outdir != '' and not os.path.isdir(outdir):
        os.makedirs(outdir)
    bad_primers = filter_primer_info(args.primer_info, args.prefix,
                                     args.Tm_min, args.Tm_max, args.Na,
                                     args.Tris, args.Mg, args.dNTPs)
            
    primers_degen = edit_primers_fasta(args.primer_degen_fasta, args.prefix,
                                       filetype='degen', to_filter=bad_primers)
    primers_expand = edit_primers_fasta(args.primer_expand_fasta, args.prefix,
                                        filetype='expand', to_filter=bad_primers)
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
