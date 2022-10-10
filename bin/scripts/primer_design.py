#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import math
import argparse
import logging
import itertools
import collections
from pprint import pprint
from subprocess import Popen, PIPE
try:
    from functools import reduce
except ImportError:
    pass
# 3rd party
from Bio.Seq import Seq
from Bio import motifs
from Bio import AlignIO
from Bio.Align import AlignInfo
import primer3

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

# argparse
desc = 'Designing primers from a set of sequences'
epi = """DESCRIPTION:
Create degenerate primers from a alignment of sequences.
primer3 is used to generate primers for a majority consensus of the alignment.
The raw primers are converted to generate primers that incorporate all sequence variation at
that alignment position.
The degenerate primers are assessed for quality, and non-optimal primers are filtered.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('target_fasta', metavar='target_fasta', type=str,
                    help='Fasta of aligned target sequences')
pg1 = parser.add_argument_group('Input/Output')
pg1.add_argument('--prefix', type=str, default='primers',
                    help='Output file prefix')
pg1.add_argument('--nontargets', type=str, default=None,
                    help='Fasta of unaligned non-target seqs')
pg2 = parser.add_argument_group('Which oligos')
pg2.add_argument('--num-raw-primers', type=int, default=100,
                 help='Number of raw primer sets to generate by primer3')
pg2.add_argument('--num-final-primers', type=int, default=10,
                 help='Max number of final primer sets.' + \
                      ' Primers with min No. of (3\') degeneracies get preference.')
pg2.add_argument('--make-oligo', type=int, default=0, choices = [0,1],
                 help='Create internal oligo?')
pg3 = parser.add_argument_group('Oligo size')
pg3.add_argument('--opt-size', type=int, default=24,
                 help='Optimum primer size')
pg3.add_argument('--min-size', type=int, default=20,
                 help='Minimum primer size')
pg3.add_argument('--max-size', type=int, default=28,
                 help='Maximum primer size')
pg4 = parser.add_argument_group('Oligo Tm')
pg4.add_argument('--opt-tm', type=float, default=60,
                 help='Optimum melting temperature')
pg4.add_argument('--min-tm', type=float, default=55,
                 help='Minimum melting temperature')
pg4.add_argument('--max-tm', type=float, default=70,
                 help='Maximum melting temperature')
pg4.add_argument('--max-tm-diff', type=float, default=1,
                 help='Maximum difference in melting temperatures between oligos')
pg4.add_argument('--oligo-DNA', type=float, default=50,
                 help='Oligo DNA concentration in the PCR [mM]')
pg4.add_argument('--dNTPs', type=float, default=0.2,
                 help='nDTPs in the PCR [mM]')
pg4.add_argument('--salt-monovalent', type=float, default=50,
                 help='Monovalent salts in the PCR [mM]')
pg4.add_argument('--salt-divalent', type=float, default=1.5,
                 help='Divalent salts in the PCR [mM]')
pg5 = parser.add_argument_group('Oligo G+C')
pg5.add_argument('--opt-gc', type=float, default=50,
                 help='Optimum G+C content')
pg5.add_argument('--min-gc', type=float, default=25,
                 help='Minimum G+C content')
pg5.add_argument('--max-gc', type=float, default=75,
                 help='Maximum G+C content')
pg5 = parser.add_argument_group('Product (amplicon) size')
pg5.add_argument('--opt-prod-size', type=int, default=250,
                    help='Optimum product (amplicon) size')
pg5.add_argument('--min-prod-size', type=int, default=150,
                    help='Minimum product (amplicon) size')
pg5.add_argument('--max-prod-size', type=int, default=350,
                    help='Maximum product (amplicon) size')
pg_int = parser.add_argument_group('Internal oligo')
pg_int.add_argument('--int-opt-size', type=int, default=20,
                    help='Internal oligo: optimal size')
pg_int.add_argument('--int-min-size', type=int, default=18,
                    help='Internal oligo: min size')
pg_int.add_argument('--int-max-size', type=int, default=22,
                    help='Internal oligo: max size')
pg_int.add_argument('--int-opt-tm', type=int, default=60,
                    help='Internal oligo: optimal size')
pg_int.add_argument('--int-min-tm', type=int, default=55,
                    help='Internal oligo: min size')
pg_int.add_argument('--int-max-tm', type=int, default=70,
                    help='Internal oligo: max size')
pg_degen = parser.add_argument_group('Degeneracy')
pg_degen.add_argument('--max-degeneracy', type=float, default=128,
                     help='Maximum degeneracy for any oligo')
pg_degen.add_argument('--max-degeneracy-3prime', type=float, default=2,
                     help='Maximum degeneracy at the 3-prime end of the fwd or rev primer')
pg_degen.add_argument('--window-3prime', type=float, default=5,
                     help='3-prime window length to check for 3-prime degeneracies')
pg_misc = parser.add_argument_group('Misc')
pg_misc.add_argument('--consensus-threshold', type=float, default=0.51,
                     help='Threshold for consensus seq construction.' + \
                     ' "N" used at all positions with no character passing the threshold.')
pg_misc.add_argument('--version', action='version', version='0.0.1')


# global variables
## IUPAC degenerate characters 
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
## reverse of IUPAC
IUPAC_R = {    
    ('A',) : 'A',
    ('C',) : 'C',
    ('G',) : 'G',
    ('T',) : 'T',
    ('A','G') : 'R',
    ('C','T') : 'Y',
    ('C','G') : 'S',
    ('A','T') : 'W',
    ('G','T') : 'K',
    ('A','C') : 'M',
    ('C','G','T') : 'B',
    ('A','G','T') : 'D',
    ('A','C','T') : 'H',
    ('A','C','G') : 'V',
    ('A','C','G','T') : 'N'
}
## The primer3 categories for oligos 
CAT_IDX = {'PRIMER_LEFT' : 'f', 'PRIMER_RIGHT' : 'r', 'PRIMER_INTERNAL' : 'i'}

# functions
def calc_consensus(aln, threshold=0.51):
    # creating consensus
    aln_sum = AlignInfo.SummaryInfo(aln)
    con_seq = aln_sum.gap_consensus(threshold=threshold, ambiguous='N')
    con_seq = str(con_seq).replace('.', 'N')
    con_seq = con_seq.replace('-', 'N')
    # check length
    con_seq_len = len(con_seq)
    aln_len = aln.get_alignment_length()
    if con_seq_len != aln_len:
        msg = 'Consensus sequence & alignment length not the same ({} != {})'
        raise ValueError(msg.format(con_seq_len, aln_len))    
    # ret
    return con_seq.lower()

def reformat_primer3_results(primers):
    """
    return: {primer_number : {primer_category : {primer3_param : value}}}
    """
    regex1 = re.compile(r'^(PRIMER_INTERNAL|PRIMER_LEFT|PRIMER_RIGHT|PRIMER_PAIR)_([0-9]+)_(.+)')
    regex2 = re.compile(r'^(PRIMER_INTERNAL|PRIMER_LEFT|PRIMER_RIGHT|PRIMER_PAIR)_([0-9]+)')
    idx = collections.defaultdict(dict)
    for k,v in primers.items():
        if k.endswith('EXPLAIN') or k.endswith('RETURNED'):
            continue
        try:
            cat,num,param = regex1.search(k).groups()
        except AttributeError:
            cat,num = regex2.search(k).groups()
            param = 'LOCATION'
        num = str(int(num)+1)
        try:
            idx[num][cat][param] = v
        except KeyError:
            idx[num][cat] = {param :  v}
    # adding internal if not present
    for num in idx.keys():
        try:
            _ = idx[num]['PRIMER_INTERNAL']
        except KeyError:
            idx[num]['PRIMER_INTERNAL'] = {
                'GC_PERCENT': None,
                'HAIRPIN_TH': None,
                'LOCATION': None,
                'PENALTY': None,
                'SELF_ANY_TH': None,
                'SELF_END_TH': None,
                'SEQUENCE': None,
                'TM': None
                }
    return idx

def calc_primers(seq, seqid = 'primers', misprime_lib=None, num_raw_primers=5,
                 opt_size=24, min_size=22, max_size=28,
                 opt_tm=60, min_tm=55, max_tm=70,
                 oligo_DNA=50, dNTPs=0.2, salt_monovalent=50, salt_divalent=1.5,
                 opt_gc=50, min_gc=25, max_gc=75,
                 opt_prod_size=300, min_prod_size=150, max_prod_size=500,
                 make_oligo = 1, int_opt_size=20, int_min_size=18, int_max_size=22,
                 int_opt_tm=60, int_min_tm=55, int_max_tm=70):
    """
    Calculating primers to with primer3-py. 
    See the primer3-py API for info on the parameters.
    Return: {primer_number : {primer_category : {primer3_param : value}}}
    """
    # dealing with short template
    max_prod_size = len(seq) if len(seq) < max_prod_size else max_prod_size
    if max_prod_size <= min_prod_size:
        msg = 'Template is too small for specified min-max product size.'
        logging.warning(msg + ' Skipping primer design')
        return None
    
    # main parameters
    params1 = {
            'SEQUENCE_ID' : seqid,
            'SEQUENCE_TEMPLATE' : seq,
            'SEQUENCE_INCLUDED_REGION' : [0,len(seq)-1]
        }
    # secondary parameters
    params2 = {
        'PRIMER_NUM_RETURN' : num_raw_primers,
        'PRIMER_PICK_INTERNAL_OLIGO' : make_oligo,
        'PRIMER_OPT_SIZE' : opt_size,
        'PRIMER_MIN_SIZE' : min_size,
        'PRIMER_MAX_SIZE' : max_size,
        'PRIMER_OPT_TM' : opt_tm,
        'PRIMER_MIN_TM' : min_tm,
        'PRIMER_MAX_TM' : max_tm,
        'PRIMER_DNA_CONC' : oligo_DNA,
        'PRIMER_DNTP_CONC' : dNTPs,
        'PRIMER_SALT_MONOVALENT' : salt_monovalent,
        'PRIMER_SALT_DIVALENT' : salt_divalent,        
        'PRIMER_OPT_GC_PERCENT' : opt_gc,
        'PRIMER_MIN_GC' : min_gc,
        'PRIMER_MAX_GC' : max_gc,
        'PRIMER_PRODUCT_OPT_SIZE' : opt_prod_size,
        'PRIMER_PRODUCT_SIZE_RANGE' : [[min_prod_size,max_prod_size]],
        'PRIMER_MAX_POLY_X' : 100,
        'PRIMER_SALT_MONOVALENT' : 50.0,
        'PRIMER_DNA_CONC' : 50.0,
        'PRIMER_MAX_NS_ACCEPTED' : 0,
        'PRIMER_MAX_SELF_ANY' : 12,
        'PRIMER_MAX_SELF_END' : 8,
        'PRIMER_PAIR_MAX_COMPL_ANY' : 12,
        'PRIMER_PAIR_MAX_COMPL_END' : 8,
        'PRIMER_INTERNAL_MAX_SELF_END' : 8,
        'PRIMER_INTERNAL_MAX_POLY_X' : 100,
        'PRIMER_INTERNAL_OPT_SIZE' : int_opt_size,
        'PRIMER_INTERNAL_MIN_SIZE' : int_min_size,
        'PRIMER_INTERNAL_MAX_SIZE' : int_max_size,
        'PRIMER_INTERNAL_OPT_TM' : int_opt_tm,
        'PRIMER_INTERNAL_MIN_TM' : int_min_tm,
        'PRIMER_INTERNAL_MAX_TM' : int_max_tm
    }
    # designing the primers
    logging.info('Designing primers...')
    primers = primer3.bindings.designPrimers(params1, params2,
                                             misprime_lib=misprime_lib)
    # re-formatting the primer3 output 
    primers = reformat_primer3_results(primers)
    msg = '  No. of raw degen primer sets: {}'
    logging.info(msg.format(len(primers.keys())))
    return primers
                
def align_loc(primers, num, cat, aln):
    """
    Determine the alignment location for the oligo
    """
    if primers[num][cat]['LOCATION'] is None:
        return None,None
    loc = list(primers[num][cat]['LOCATION'])
    aln_len = aln.get_alignment_length()
    if cat == 'PRIMER_RIGHT':
        start = loc[0] - loc[1] + 1
        end = loc[0] + 1
    else:
        start = loc[0]
        end = loc[0] + loc[1]
    return start,end

def get_subseqs(aln, start, end, cat, allow_gaps=False):
    """
    Parsing out sub-sequences from a Bio.Align object.
    """
    subseqs = []
    for aln_seq in aln:
        seq = aln_seq[start:end].seq
        if cat == 'PRIMER_RIGHT':
            seq = seq.reverse_complement()
        if allow_gaps is False and str(seq).find('-') != -1:
            return None
        subseqs.append(seq.upper())
    return subseqs
    
def seqs2degen(seqs):
    """
    Convert a list of sequences into a degenerate consensus sequences.
    Asumming that all of the sequences are aligned (assessing position by position).
    Returning upper case degen sequence.
    """
    if seqs is None:
        return ''
    # converting to list of lists (individual chars)
    seqs = [[y for y in str(x).upper()] for x in seqs]
    seq_lens = set([len(x) for x in seqs])
    if not len(seq_lens) == 1:
        raise KeyError('Sequences are not the same length!')
    # getting seq variation at each position, then mapping to IUPAC of degen chars
    degen_seq = ''
    for i in range(list(seq_lens)[0]):
        chars = set([seq[i] for seq in seqs])
        if '-' in chars:
            return ''
        if 'N' in chars:
            chars = {'A', 'C', 'G', 'T'}
        try:
            degen_seq += IUPAC_R[tuple(sorted(set(chars)))]
        except KeyError:
            msg = 'Cannot map chars to degen character: {}'
            raise KeyError(msg.format(','.join(list(chars))))
    # checking the degen seq_len = input_seqs
    if len(degen_seq) != list(seq_lens)[0]:
        raise ValueError('Degenerate seq len != non-degen seqs!')
    return degen_seq

def expand_degen(seq, degen_max, window=None):
    """
    From string of nucleotides, create list of sequences 
    comprising all non-degenerate versions of the original sequence.
    """
    seq = [IUPAC[bp] for bp in str(seq)]
    if window is not None:
        window = -window if window > 0 else window
        seq = seq[int(window):]
    degen = reduce((lambda x, y: x * y), [len(x) for x in seq])
    if degen > degen_max:
        return None     
    return [''.join(x) for x in itertools.product(*seq)]

def calc_amplicon(primers, aln, seq_cats):
    """
    Calculating the average amplicon size by:
      - pulling out subseqs of the alignment
      - removing any gaps in the subseqs
      - calculating length of the subseqs
    """
    for num in primers.keys():
        # getting primer start-end info
        amplicon_lens = []
        left_start,left_end = [],[]
        right_start,right_end = [],[]
        for cat in ['PRIMER_LEFT','PRIMER_RIGHT']:
            for seq in primers[num][cat].keys():
                # getting start-end
                if cat == 'PRIMER_LEFT':
                    left_start.append(primers[num][cat][seq]['start'])
                    left_end.append(primers[num][cat][seq]['end'])
                elif cat == 'PRIMER_RIGHT':
                    right_start.append(primers[num][cat][seq]['start'])
                    right_end.append(primers[num][cat][seq]['end'])
                else:
                    raise ValueError
        # determining amplicon seq & length
        for i in range(len(left_start)):
            if left_end[i] <= right_start[i]:
                start,end = left_start[i],right_end[i]
            else:
                start,end = right_start[i],left_end[i]
            subseqs = get_subseqs(aln, start, end, 'PRIMER_LEFT', allow_gaps=True)
            for seq in subseqs:
                seq = str(seq).replace('-', '').replace('.', '')
                amplicon_lens.append(len(seq))
        # averaging lengths
        amp_len_avg = avg(amplicon_lens)
        amp_len_sd = sd(amplicon_lens)
        # adding amplicon info to primer dict
        for cat in seq_cats:
            for seq in primers[num][cat].keys():
                primers[num][cat][seq]['amplicon_size_avg'] = amp_len_avg
                primers[num][cat][seq]['amplicon_size_sd'] = amp_len_sd
    return primers

def filter_incomplete_primer_sets(primers, seq_cats):
    """
    Filtering incomplete primer sets
    """
    # any primer sets lacking a seq_cat are filtered
    missing_cnt = {x:0 for x in seq_cats}
    to_rm = set()
    for num in primers.keys():
        for cat in seq_cats:
            try:
                _ = primers[num][cat]
            except KeyError:
                to_rm.add(num)
                missing_cnt[cat] += 1
    ## filtering
    for x in to_rm:
        primers.pop(x, None)
    # status
    msg = '  No. of incomplete primer sets due to oligo filtering: {}'
    logging.info(msg.format(len(to_rm)))
    for k in seq_cats:
        msg = '    No. of {} oligo missing: {}'
        logging.info(msg.format(k, missing_cnt[k]))
    msg = '  No. of primers retained: {}'
    logging.info(msg.format(len(primers.keys())))
    return primers

def get_len(x):
    try:
        return len(x)
    except TypeError:
        return 0

def calc_degen(primers, aln, degen_max, degen_3prime_max, window_3prime=5,
               internal_oligo=0):
    """
    Calculating degenerate primers
    """
    logging.info('Calculating degenerate primers...')
    if primers is None:
        return None
    seq_cats = {'PRIMER_LEFT', 'PRIMER_RIGHT'}
    if internal_oligo == 1:
        seq_cats.add('PRIMER_INTERNAL')
    status = {'no location' : 0,
              'no primary subseqs' : 0,
              'oligo degen exceeded' : 0,
              '3p max degen exceeded' : 0} 
    degen_primers = collections.defaultdict(dict)
    for num in primers.keys():
        for cat in primers[num].keys():
            if not cat in seq_cats:
                continue
            # location in the alignemnt, depending on the olig
            if primers[num][cat]['LOCATION'] is None:
                status['no location'] += 1
                continue
            start,end = align_loc(primers, num, cat, aln)
            # subsetting target sequences to just primer overlap
            subseqs = get_subseqs(aln, start, end, cat)
            if subseqs is None:
                status['no primary subseqs'] += 1
                continue
            # calculating degnerate sequence across all targets
            degen_seq = seqs2degen(subseqs)            
            # checking for 3p max degen
            if cat.endswith('LEFT') or cat.endswith('RIGHT'):
                degen_3p = expand_degen(degen_seq, degen_3prime_max,
                                        window=window_3prime)
                if degen_3p is None:
                    status['3p max degen exceeded'] += 1
                    continue
            # expanded sequences based on degeneracies 
            expand_seq = expand_degen(degen_seq, degen_max)
            if expand_seq is None:
                status['oligo degen exceeded'] += 1
                continue
            # summing up info
            stats = {
                'length' : len(degen_seq),
                'expanded' : expand_seq,
                'degeneracy' : get_len(expand_seq),
                'degeneracy_3p' : get_len(degen_3p),
                'start' : start,
                'end' : end,
                'amplicon_size_consensus' : primers[num]['PRIMER_PAIR']['PRODUCT_SIZE']
                }
            try:                
                degen_primers[num][cat][str(degen_seq)] = stats
            except KeyError:
                degen_primers[num][cat] = {str(degen_seq) : stats}
    # status
    for k,v in status.items():
        logging.info('  Filtered {} oligos due to {}'.format(v,k))
    msg = '  No. of (partial) primer sets remaining: {}'
    logging.info(msg.format(len(degen_primers.keys())))
    # filter out partial primer sets
    degen_primers = filter_incomplete_primer_sets(degen_primers, seq_cats)
    # calculating amplicon size based on alignment
    degen_primers = calc_amplicon(degen_primers, aln, seq_cats)
    # return
    return degen_primers

def avg(x):
    if len(x) == 0:
        return 0
    else:
        return sum(x) / float(len(x))

def sd(x):
    if len(x) == 0:
        return 0
    mu = avg(x)
    sigma_sq = sum([(y - mu)**2 for y in x]) / float(len(x))
    return math.sqrt(sigma_sq)

def calc_GC(seq):
    """
    Calc. seq G+C content
    """
    gc = {'G', 'C', 'g', 'c'}
    return len([x for x in seq if x in gc]) / float(len(seq)) * 100

def expanded_primer_stats(degen, oligo_DNA=50, dNTPs=0.2,
                          salt_monovalent=50, salt_divalent=1.5):
    """
    Calculating per-non-degen-primer stats (Tm) and averaging
    """
    logging.info('Calculating stats on primer sets...')
    if degen is None:
        return None
    for num in degen.keys():
        for cat in degen[num].keys():
            for degen_seq in degen[num][cat].keys():
                stats = {'Tm' : [],
                         'GC' : [],
                         'hairpin' : [],
                         'homodimer' : []}
                # stats on each expanded primer
                for seq in list(degen[num][cat][degen_seq]['expanded']):
                    # degeneracies
                    # melting temp
                    stats['Tm'].append(primer3.calcTm(seq,
                                                      dna_conc=oligo_DNA,
                                                      dntp_conc=dNTPs,
                                                      mv_conc=salt_monovalent,
                                                      dv_conc=salt_divalent))
                    # GC
                    stats['GC'].append(calc_GC(seq))
                    # hairpin
                    stats['hairpin'].append(primer3.calcHairpin(seq,
                                                                dna_conc=oligo_DNA,
                                                                dntp_conc=dNTPs,
                                                                mv_conc=salt_monovalent,
                                                                dv_conc=salt_divalent).tm)
                    # homodimer
                    stats['homodimer'].append(primer3.calcHomodimer(seq,
                                                                    dna_conc=oligo_DNA,
                                                                    dntp_conc=dNTPs,
                                                                    mv_conc=salt_monovalent,
                                                                    dv_conc=salt_divalent).tm)
                # summarizing stats (average & std)
                for k,v in stats.items():
                    degen[num][cat][degen_seq][k] = [avg(v), sd(v)]
                                                            
    return degen    

def filter_primer_sets_by_tm_diff(primers, max_tm_diff=1):
    max_tm_diff = float(max_tm_diff)
    to_rm = set()
    for num in primers.keys():
        # getting Tm values for all oligos in primer set
        Tm_list = []
        for cat in primers[num].keys():
            for seq in primers[num][cat].keys():
                Tm_list.append(primers[num][cat][seq]['Tm'][0])
        # checking min/max diff
        if max(Tm_list) - min(Tm_list) > max_tm_diff:
            to_rm.add(num)
    ## filtering
    for x in to_rm:
        primers.pop(x, None)
    # status
    logging.info('  No. primer sets exceeding max-tm-diff: {}'.format(len(to_rm)))
    logging.info('  No. of retained primers: {}'.format(len(primers.keys())))
    return primers

def subset_primer_sets(primers, max_final_primers=100):
    """
    primers: {primer_set_id : oligo_cat : {seq : {characteristic : value}}}
    """
    max_final_primers = int(max_final_primers)
    if len(primers.keys()) <= max_final_primers:
        return primers
    else:
        n_to_filter = len(primers.keys()) - max_final_primers
    msg = 'Filtering {} primer sets to {} final sets...'
    logging.info(msg.format(len(primers.keys()), max_final_primers))

    # getting degen of each primer set
    degen_stats = {}
    for num in primers.keys():
        # getting Tm values for all oligos in primer set
        degen_all = []
        degen_3prime = []
        for cat in primers[num].keys():
            for seq in primers[num][cat].keys():
                degen_all.append(primers[num][cat][seq]['degeneracy'])
                degen_3prime.append(primers[num][cat][seq]['degeneracy_3p'])
        degen_stats[num] = [sum(degen_all),sum(degen_3prime)]
    # filtering those with "worst" degen
    to_rm = set()    
    for num,stats in sorted(degen_stats.items(), key=lambda x: (x[1][0], x[1][1]), reverse=True):
        if len(to_rm) >= n_to_filter:
            break
        to_rm.add(num)    
    ## actual filtering
    for x in to_rm:
        primers.pop(x, None)
    # status
    logging.info('  No. of retained primers: {}'.format(len(primers.keys())))
    return primers

def write_primer_fasta(primers, prefix, filetype='degen'):
    """
    Write primer sequences in fasta format.
    'degen' = degenerate primers written
    'expand' = expanded (non-degen) primers written
    """
    if primers is None:
        primers = {}
    outfile = prefix + '_{}.fna'.format(filetype)
    with open(outfile, 'w') as outF:
        for num in primers.keys():
            for cat in primers[num].keys():
                fwi = CAT_IDX[cat]
                for degen_seq in primers[num][cat].keys():
                    if filetype == 'degen':
                        outF.write('>{}{}\n'.format(num, fwi))
                        outF.write(degen_seq + '\n')
                    elif filetype == 'expand':
                        for i,exp_seq in enumerate(primers[num][cat][degen_seq]['expanded']):
                            outF.write('>{}{}_{}\n'.format(num, fwi, i))
                            outF.write(exp_seq + '\n')                            
                    else:
                        raise ValueError                        
    logging.info('File written: {}'.format(outfile))

def write_primer_info(primers, prefix):
    """
    Write tsv of primer metadata
    """
    if primers is None:
        primers = {}
    outfile = prefix + '.tsv'
    header = ['primer_set', 'amplicon_size_consensus',
              'amplicon_size_avg', 'amplicon_size_sd', 
              'primer_id', 'primer_type', 'sequence',      
              'length', 'degeneracy', 'degeneracy_3prime',
              'position_start', 'position_end',
              'Tm_avg', 'Tm_sd', 'GC_avg', 'GC_sd',
              'hairpin_avg', 'hairpin_sd',
              'homodimer_avg', 'homodimer_sd']
    with open(outfile, 'w') as outF:
        outF.write('\t'.join(header) + '\n')
        for num in sorted(primers.keys(), key=lambda x: int(x)):
            for cat in primers[num].keys():
                fwi = CAT_IDX[cat]
                for degen_seq in primers[num][cat].keys():
                    x = [
                        str(num),
                        str(primers[num][cat][degen_seq]['amplicon_size_consensus']),
                        str(primers[num][cat][degen_seq]['amplicon_size_avg']),
                        str(primers[num][cat][degen_seq]['amplicon_size_sd']),
                        '{}{}'.format(num, fwi),
                        cat,
                        degen_seq,
                        str(primers[num][cat][degen_seq]['length']),
                        str(primers[num][cat][degen_seq]['degeneracy']),
                        str(primers[num][cat][degen_seq]['degeneracy_3p']),
                        str(primers[num][cat][degen_seq]['start']),
                        str(primers[num][cat][degen_seq]['end']),
                        str(primers[num][cat][degen_seq]['Tm'][0]),
                        str(primers[num][cat][degen_seq]['Tm'][1]),
                        str(primers[num][cat][degen_seq]['GC'][0]),
                        str(primers[num][cat][degen_seq]['GC'][1]),
                        str(primers[num][cat][degen_seq]['hairpin'][0]),
                        str(primers[num][cat][degen_seq]['hairpin'][1]),
                        str(primers[num][cat][degen_seq]['homodimer'][0]),
                        str(primers[num][cat][degen_seq]['homodimer'][1]),  
                    ]
                    outF.write('\t'.join(x) + '\n')
    logging.info('File written: {}'.format(outfile))
    
def main(args):
    """
    Main interface
    """
    # output directory
    outdir = os.path.split(args.prefix)[0]
    if outdir != '' and not os.path.isdir(outdir):
        os.makedirs(outdir)
    # calculate consensus for targets
    target_aln = AlignIO.read(args.target_fasta, 'fasta')
    con_seq = calc_consensus(target_aln, threshold=args.consensus_threshold)
    # calculate primers
    primers = calc_primers(con_seq,
                           num_raw_primers=args.num_raw_primers,
                           opt_size=args.opt_size,
                           min_size=args.min_size,
                           max_size=args.max_size,
                           opt_tm=args.opt_tm,
                           min_tm=args.min_tm,
                           max_tm=args.max_tm,
                           oligo_DNA=args.oligo_DNA,
                           dNTPs=args.dNTPs,
                           salt_monovalent=args.salt_monovalent,
                           salt_divalent=args.salt_divalent,
                           opt_gc=args.opt_gc,
                           min_gc=args.min_gc,
                           max_gc=args.max_gc,
                           opt_prod_size=args.opt_prod_size,
                           min_prod_size=args.min_prod_size,
                           max_prod_size=args.max_prod_size,
                           make_oligo=args.make_oligo,
                           int_opt_size=args.int_opt_size,
                           int_min_size=args.int_min_size,
                           int_max_size=args.int_max_size,
                           int_opt_tm=args.int_opt_tm,
                           int_min_tm=args.int_min_tm,
                           int_max_tm=args.int_max_tm)
    
    # generating degenerate primers from the cons-seq primers
    degen = calc_degen(primers, target_aln, args.max_degeneracy,
                       degen_3prime_max=args.max_degeneracy_3prime,
                       window_3prime=args.window_3prime,
                       internal_oligo=args.make_oligo)
    # calculating parameters for "expanded" primers
    degen = expanded_primer_stats(degen,
                                  oligo_DNA=args.oligo_DNA,
                                  dNTPs=args.dNTPs,
                                  salt_monovalent=args.salt_monovalent,
                                  salt_divalent=args.salt_divalent)
    # filtering based on oligo Tm differences
    degen = filter_primer_sets_by_tm_diff(degen, max_tm_diff=args.max_tm_diff)
    # subset to max number of primers
    degen = subset_primer_sets(degen, max_final_primers=args.num_final_primers)    
    # write primer fasta
    ## degenerate
    write_primer_fasta(degen, args.prefix, filetype='degen')
    ## expanded (all non-degenerate)
    write_primer_fasta(degen, args.prefix, filetype='expand')
    # write primer info
    write_primer_info(degen, args.prefix)
    
    

if __name__ == '__main__':
    args = parser.parse_args()    
    main(args)
