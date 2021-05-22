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

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

# arg parse
desc = 'Filter out undesireable primers'
epi = """DESCRIPTION:
Filtering primerprospector::generate_primers_denovo.py output.
Primers can be filtered on certain criteria.
The fitlered primers are then sorted & condensed to degenerate primers.
The degenerate primers are filtered by --degen-max and expanded to all non-degen sequences.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('primer_file', metavar='primer_file', type=str,
                    help='primerprospector::generate_primers_denovo.py output file(s)')
parser.add_argument('--prefix', type=str, default='primers_filtered',
                    help='Output file prefix')
parser.add_argument('--n-primers-max', type=int, default=0,
                    help='Max number of degenerate primers kept (0 = all kept)')
parser.add_argument('--amp-len-min', type=int, default=200,
                    help='Min amplicon length of the primer pair')
parser.add_argument('--amp-len-max', type=int, default=800,
                    help='Max amplicon length of the primer pair')
parser.add_argument('--degen-max', type=int, default=64,
                    help='Max no. of degeneracies in the primer')
parser.add_argument('--gc-min', type=float, default=40,
                    help='Min GC of the primer pair')
parser.add_argument('--gc-max', type=float, default=60,
                    help='Max GC of the primer pair')
parser.add_argument('--tmp-dir', type=str, default='primer_filter_tmp',
                    help='Temporary file directory')
parser.add_argument('--tmp-dir-keep', action='store_true', default=False,
                    help='Keep temporary file directory')
parser.add_argument('--version', action='version', version='0.0.1')


def load_primers(infile):
    """
    Reading in primerprospector output
    """
    logging.info('Reading file: {}'.format(infile))
    primer_cnt = 0
    primers = collections.defaultdict(dict)
    with open(infile) as inF:
        for line in inF:
            if line.startswith('#'):
                continue
            line = line.rstrip().split(',')
            if line[0] == 'C':
                primer_cnt += 1
                primers[primer_cnt]['stats'] = line
            else:
                try:
                    primers[primer_cnt]['primers'].append(line)
                except KeyError:
                    primers[primer_cnt]['primers'] = [line]
    logging.info('  No. of raw primers: {}'.format(len(primers.keys())))
    return primers

def calc_perc_gc(seq):
    GC = ('G', 'C')
    GC = len([x for x in seq if x in GC])
    return float(GC) / len(seq) * 100

def filter_primers(primers, gc_min=40, gc_max=60, amp_len_min=200, amp_len_max=800):
    """
    Filtering primer dict
    """
    logging.info('Initial filtering of primers...')
    to_rm = []
    status = {'gaps' : 0, 'gc' : 0, 'amp_len' : 0}
    for primer_id,d in primers.items():
        to_filter = False        
        for i,primer in enumerate(d['primers']):
            # filtering gaps
            if '-' in primer[2] or '-' in primer[3]:
                to_filter = True
                status['gaps'] += 1
                break
            # amplicon
            if int(primer[1]) < amp_len_min or int(primer[1]) > amp_len_max:
                to_filter = True
                status['amp_len'] += 1
                break
            # filtering GC content
            if i == 0:
                GC = calc_perc_gc(primer[2] + primer[3])
                if GC < gc_min or GC > gc_max:
                    to_filter = True
                    status['gc'] += 1
                    break
        if to_filter:
            to_rm.append(primer_id)
    for x in to_rm:
        primers.pop(x, None)
    logging.info('  No. of primers filtered: {}'.format(len(to_rm)))
    logging.info('    Contained gaps: {}'.format(status['gaps']))
    logging.info('    G+C out of bounds: {}'.format(status['gc']))
    logging.info('    Amplicon len. out of bounds: {}'.format(status['amp_len']))
    logging.info('  No. of remaining primers: {}'.format(len(primers.keys())))    

def write_primer_csv(primers, outdir):
    if outdir != '' and not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, 'filtered_primers.csv')
    with open(outfile, 'w') as outF:
        for primer_id,d in primers.items():
            outF.write(','.join(d['stats']) + '\n')
            for primer in d['primers']:
                outF.write(','.join(primer) + '\n')
    logging.info('File written: {}'.format(outfile))
    return outfile

def sort_denovo_primers(exe, primer_file, tmp_dir):
    """
    Calling sort_denovo_primers.py
    """
    logging.info('Sorting primers...')
    cmd = [exe, '-i', primer_file, '-o', tmp_dir]
    logging.info('  CMD: {}'.format(' '.join(cmd)))
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    rc = p.returncode
    if rc != 0:
        for line in output.decode().split('\n'):
            sys.stderr.write(line + '\n')
        raise ValueError(err)
    files = [os.path.join(tmp_dir, 'formatted_primers.txt'),
             os.path.join(tmp_dir, 'primers_details.txt')]
    for F in files:
        logging.info('  File written: {}'.format(F))
    return files

IUPAC = {
    'A' : ('A'),
    'C' : ('C'),
    'G' : ('G'),
    'T' : ('T'),
    'R' : ('A','G'),
    'Y' : ('C','T'),
    'S' : ('G','C'),
    'W' : ('A','T'),
    'K' : ('G','T'),
    'M' : ('A','C'),
    'B' : ('C','G','T'),
    'D' : ('A','G','T'),
    'H' : ('A','C','T'),
    'V' : ('A','C','G'),
    'N' : ('A','C','G','T')
    }

def _exp_degen(seq, degen_max):
    """
    From string of nucleotides, create list of sequences 
    comprising all non-degenerate versions of the original sequence.
    """
    seq = [IUPAC[bp] for bp in seq]
    degen = reduce((lambda x, y: x * y), [len(x) for x in seq])
    if degen > degen_max:
        return None
    return [''.join(x) for x in itertools.product(*seq)]        

def read_primers(infile, primer='f', degen_max=128, n_primers_max=10):
    """
    Reading in primers
    """
    logging.info('Reading file: {}'.format(infile))
    primers = collections.defaultdict(dict)
    with open(infile) as inF:
        for line in inF:
            if line.startswith('#'):
                continue
            line = line.rstrip().split('\t')
            if line[0] == '' or not line[0].endswith(primer):
                continue
            # expanding degen sequences
            expanded_primers = _exp_degen(line[1], degen_max)
            primers[line[0]]['degen'] = [line[1]]
            primers[line[0]]['expand'] = expanded_primers
            # max number of primer cutoff
            if n_primers_max > 0 and len(primers.keys()) >= n_primers_max:
                break
    return primers

def filter_by_degen(primers, degen_max):
    to_rm = set()
    for primer,d in primers.items():
        if d['f']['expand'] is None or d['r']['expand'] is None:
            to_rm.add(primer)
        elif len(d['f']['expand']) + len(d['r']['expand']) > degen_max:
            to_rm.add(primer)
    for x in to_rm:
        primers.pop(x, None)
    logging.info('  No. of primer pairs exceeding degen-max: {}'.format(len(to_rm)))

def expand_degenerate(infile, degen_max, n_primers_max):
    """
    Loading sorted primers and expanding degeneracies in sequences
    """
    logging.info('Loading file: {}'.format(infile))
    primers = collections.defaultdict(dict)
    for primer,d in read_primers(infile, 'f', degen_max, n_primers_max).items():
        x = primer.rstrip('f')
        primers[x]['f'] = d
    for primer,d in read_primers(infile, 'r', degen_max, n_primers_max).items():
        x = primer.rstrip('r')
        primers[x]['r'] = d
    for primer,d in primers.items():
        try:
            _ = d['f']
            _ = d['r']
        except KeyError:
            msg = 'Primer does not have fwd + rev: {}'
            raise KeyError(msg.format(primer))
    logging.info('Filtering primers by degeneracy...')
    n_filtered = filter_by_degen(primers, degen_max)
    
    logging.info('  No. of degen-allowed primer pairs: {}'.format(len(primers.keys())))
    fwd_all = sum([len(x['f']['expand']) for x in primers.values()])
    rev_all = sum([len(x['r']['expand']) for x in primers.values()])
    logging.info('  No. of degen-expanded primers: {}'.format(fwd_all + rev_all))
    return primers

def _write_primers_fasta(d, primer_id, seqtype, outF, fr):
    for i,seq in enumerate(d[fr][seqtype]):
        if seqtype == 'expand':
            seq = '>{}{}_{}\n{}\n'.format(primer_id, fr, i, seq)
        else:
            seq = '>{}{}\n{}\n'.format(primer_id, fr, seq)
        outF.write(seq)    

def write_primers_fasta(primers, prefix, seqtype='expand'):
    outfile = prefix + '_{}.fna'.format(seqtype)
    with open(outfile, 'w') as outF:
        for primer_id,d in primers.items():
            _write_primers_fasta(d, primer_id, seqtype, outF, 'f')
            _write_primers_fasta(d, primer_id, seqtype, outF, 'r')
    logging.info('File written: {}'.format(outfile))
    
def write_primers_metadata(primers, primer_info, sorted_primer_file, prefix):
    logging.info('Reading file: {}'.format(sorted_primer_file))
    header = ['Conserved_Xmer', 'Unaligned_Index', 'Aligned_Index',
              'Number_Hits', 'Percent_Match', 'Nonspecific_Match',
              '5p_Standard_Index', '3p_Standard_Index',
              'F_Primer_IUPAC_Sequence', 'F_Primer_Filtered_Sequence',
              'F_Primer_Majority_Consensus', 'R_Primer_IUPAC_Sequence',
              'R_Primer_Filtered_Sequence', 'R_Primer_Majority_Consensus']
    primer_metadata = {}   # {fwd_seq : metadata}
    with open(sorted_primer_file) as inF:
        for line in inF:
            if line.startswith('#'):
                continue
            line = line.rstrip().split(',')
            if line[0] == '':
                continue
            line = [x if x != '' else 'N/A' for x in line]
            primer_metadata[line[8]] = line
    # finding primers of interest
    outfile = prefix + '.tsv'
    with open(outfile, 'w') as outF:
        outF.write('\t'.join(['primer_id'] + header) + '\n')
        for primer,d in primers.items():
            seq = d['f']['degen'][0]
            try:
                record = primer_metadata[seq]
            except KeyError:
                msg = 'Cannot find primer in metadata: {}'
                logging.info(msg.format(seq))
            outF.write('\t'.join([primer] + record) + '\n')
    logging.info('File written: {}'.format(outfile))
    
def main(args):    
    # checking for primerprospector exe
    exe = 'sort_denovo_primers.py'
    if find_executable(exe) is None:
        msg = 'Cannot find executable: {}'
        raise ValueError(msg.format(exe))
    # reading primer file
    primer_info = load_primers(args.primer_file)
    # filtering of primers
    filter_primers(primer_info,
                   args.gc_min, args.gc_max,
                   args.amp_len_min, args.amp_len_max)
    primer_file = write_primer_csv(primer_info, args.tmp_dir)    
    # running sort_denovo_primers
    sorted_primer_files = sort_denovo_primers(exe, primer_file, args.tmp_dir)
    # loading formatted primers & expanding degenerate sequences
    primers = expand_degenerate(sorted_primer_files[0], args.degen_max, args.n_primers_max)
    # writing primers
    write_primers_fasta(primers, args.prefix, seqtype='degen')
    write_primers_fasta(primers, args.prefix, seqtype='expand')
    # writing metadata
    write_primers_metadata(primers, primer_info, sorted_primer_files[1], args.prefix)
    # clean up
    if args.tmp_dir_keep is False and os.path.isdir(args.tmp_dir):
        shutil.rmtree(args.tmp_dir)
        logging.info('tmp-dir removed: {}'.format(args.tmp_dir))        
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
