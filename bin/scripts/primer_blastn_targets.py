#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
import collections

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'blastn of primers versus the target gene sequences'
epi = """DESCRIPTION:
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('blast_hits', metavar='blast_hits', type=str,
                    help='tab-delim table of blast hits')
parser.add_argument('primer_degen_fasta', metavar='primer_degen_fasta', type=str,
                    help='Fasta file of degenerate primer sequences')
parser.add_argument('primer_expand_fasta', metavar='primer_expand_fasta', type=str,
                    help='Fasta file of expanded primer sequences')
parser.add_argument('primer_info', metavar='primer_info', type=str,
                    help='tab-delim table of primer info')
parser.add_argument('--prefix', type=str, default='blast_filtered_primers',
                    help='Output file prefix')
parser.add_argument('--perc-len', type=float, default=80,
                    help='Hit is not legit if (aln_len/query_len*100) < X')
parser.add_argument('--min-amplicon-len', type=int, default=100,
                    help='Min. amplicon size to be considered an off-target hit')
parser.add_argument('--max-amplicon-len', type=int, default=1500,
                    help='Max. amplicon size to be considered an off-target hit')
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

def to_pos_strand(sstart, send, slen):
    """
    If blast position info is scaled on the negative strand, flipping the info to the positive
    """
    # nothing required if positive strand
    if send >= sstart:
        return sstart,send,'+'
    # convert to positive strand position if on negative
    sstart = slen - sstart + 1
    send = slen - send + 1
    if sstart < 1 or send < 0:
        msg = 'flipping position to + strand resulted in value of < 1'
        raise ValueError(msg)
    return sstart,send,'-'

def load_blast_results(infile, perc_len=90.0):
    """
    Expected input format: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen
    
    Return: {primer_id : {fwd|rev : [[saccver, sstart, send],...]}}
    """
    logging.info('Reading file: {}'.format(infile))
    regex = re.compile(r'^([0-9]+)([fr])_([0-9]+)$')
    status = {'target' : 0, 'non-sig' : 0}
    blast_hits = collections.defaultdict(dict)
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            # filter by perc. alignment length
            aln_len = float(line[3]) / float(line[-1]) * 100
            if aln_len < perc_len:
                status['non-sig'] += 1
                continue
            # keeping primer
            try:
                primer_id,fr,derep_id = regex.search(line[0]).groups()
            except AttributeError:
                msg = 'Primer name formatted incorrectly: {}'
                raise AttributeError(msg.format(line[0]))
            ## adding hit to index
            saccver = line[1]
            slen = int(line[-2])
            sstart = int(line[8])
            send = int(line[9])
            sstart,send,strand = to_pos_strand(sstart, send, slen)
            try:
                blast_hits[primer_id][saccver][fr].append([sstart, send, strand])
            except KeyError:
                try:
                    blast_hits[primer_id][saccver][fr] = [[sstart, send, strand]]
                except KeyError:
                    blast_hits[primer_id][saccver] = {fr : [[sstart, send, strand]]}
            status['target'] += 1
    logging.info('  No. of signfiicant hits: {}'.format(status['target']))
    logging.info('  No. of non-significant hits: {}'.format(status['non-sig']))
    return blast_hits

def calc_amp_len(x_start, x_end, y_start, y_end):
    """
    Determining the amplicon length based on fwd-rev primer blast hits
    """
    x_min = min([x_start, x_end])
    x_max = max([x_start, x_end])
    y_min = min([y_start, y_end])
    y_max = max([y_start, y_end])    
    if y_min > x_max:
        return y_max - x_min
    else:
        return x_max - y_min

def avg(x):
    try:
        y = sum(x) / len(x)
    except ZeroDivisionError:
        y = 0
    return y

def sd(x):
    mu = avg(x)
    if mu == 0:
        return 0
    try:
        y = sum([(z - mu)**2 for z in x]) / len(x)
    except ZeroDivisionError:
        y = 0
    return y    

def hits2info(blast_hits, min_amplicon_len, max_amplicon_len):
    """
    Converting the blastn hit info to primer annealing position & amplicion length info
    Return: {primer_id : [fwd_ave_3prime_pos, rev_ave_3prime_pos, amplicon_ave_len]}
    """
    logging.info('Summarizing blast info...')
    info = {}
    status = {'no_amp' : 0}
    for primer_id in blast_hits.keys():
        amp_lens = []
        fwd_3prime = []
        rev_3prime = []
        for saccver,dd in blast_hits[primer_id].items():
            # hits for primer pair?
            try:
                _,_ = dd['f'],dd['r']
            except KeyError:
                continue
            # amplicon size of any expanded fwd-rev pair            
            for x in dd['f']:
                for y in dd['r']:
                    amp_len = calc_amp_len(x[0], x[1], y[0], y[1])
                    if (x[2] != y[2] and amp_len >= min_amplicon_len
                        and amp_len <= max_amplicon_len):
                        amp_lens.append(amp_len)
                        fwd_3prime.append(x[1] if x[2] == '+' else x[0])
                        rev_3prime.append(y[1] if y[2] == '+' else y[0])
        # per-primer summary
        f3p = avg(fwd_3prime)
        r3p = avg(rev_3prime)
        amp_strand = '+' if r3p > f3p else '-'
        info[primer_id] = [round(avg(amp_lens),1), round(sd(amp_lens),1),
                           round(f3p,1), round(r3p,1), amp_strand]
        if avg(amp_lens) == 0:
            status['no_amp'] += 1
    # summary
    logging.info('  No. of primers: {}'.format(len(info.keys())))
    logging.info('  No. of primers with no relevant amplicon: {}'.format(status['no_amp']))
    return info    

def load_primer_fasta(infile, nonspec_primers, seqtype='degen'):
    """
    Loading primer fasta. Excluding non-specific primers (nonspec_primers).    
    """
    logging.info('Reading file: {}'.format(infile))
    seqs = collections.defaultdict(dict)
    regex = re.compile(r'^[0-9]+([rf])')
    seqid = None
    n_seq = 0
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            if line.startswith('>'):
                seqid = line.lstrip('>')
                fr = regex.search(seqid).groups()[0]
                n_seq += 1
                continue
            # primer id 
            primer_id = seqid.split('_')[0].rstrip(fr)
            if seqtype == 'degen':
                num = 0
            elif seqtype == 'expand':
                num = int(seqid.split('_')[1])
            else:
                raise ValueError('seqtype not recognized')
            # adding sequences
            try:
                seqs[primer_id][fr][num] += line
            except KeyError:
                try:
                    seqs[primer_id][fr][num] = line
                except KeyError:
                    seqs[primer_id][fr] = {num : line}
    logging.info('  No. of total primer pairs: {}'.format(len(seqs.keys())))
    logging.info('  No. of total primer seqs: {}'.format(n_seq))
    # filtering
    for p in nonspec_primers:
        seqs.pop(p, None)
    logging.info('  No. of retained primers: {}'.format(len(seqs.keys())))
    return seqs

def write_primer_fasta(primers, prefix, seqtype='degen'):
    """
    Writing primer fasta file
    """
    outdir = os.path.split(prefix)[0]
    if outdir != '' and not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfile = prefix + '_{}.fna'.format(seqtype)
    cnt = 0
    with open(outfile, 'w') as outF:
        for primer,d in primers.items():
            for fr,dd in d.items():
                for num in sorted(dd.keys()):
                    if seqtype == 'degen':
                        outF.write('>{}{}\n{}\n'.format(primer, fr, d[fr][num]))
                    elif seqtype == 'expand':
                        outF.write('>{}{}_{}\n{}\n'.format(primer, fr, num, d[fr][num]))
                    else:
                        raise ValueError('seqtype not recognized')
                    cnt += 1
    logging.info('File written: {}'.format(outfile))
    logging.info('  No of primer pairs written: {}'.format(cnt))

def edit_primer_info(infile, new_info, primers, prefix):
    """
    Loading info on primers, editing, and writing
    """
    if infile is None:
        return Nnoe
    outdir = os.path.split(prefix)[0]
    if outdir != '' and not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfile = prefix + '.tsv'
    
    logging.info('Reading file: {}'.format(infile))
    header = {}
    cnt = 0
    with open(infile) as inF, open(outfile, 'w') as outF:
        for i,line in enumerate(inF):
            line = line.rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                to_add = ['ave_amplicon_len', 'sd_amplicon_len',
                          'fwd_3prime_position_plus_strand',
                          'rev_3prime_position_plus_strand',
                          'amplicon_strand']
                outF.write('\t'.join(line + to_add) + '\n')
                continue
            primer_id = str(line[header['primer_id']])
            # filtering based on provided list of primers
            try:
                x = primers[primer_id]
            except KeyError:
                continue
            # getting new primer info
            try:
                y = new_info[primer_id]
            except KeyError:
                msg = 'Cannot find "{}" in new primer info'
                raise KeyError(msg.format(primer_id))
            if len(x.keys()) > 0:
                outF.write('\t'.join(line + [str(z) for z in y]) + '\n')
                cnt += 1
    logging.info('  File written: {}'.format(outfile))            
    logging.info('  No. of records written: {}'.format(cnt))
    
def main(args):
    """
    * parse blastn results to get location of primer annealing & amplicon length
    """
    # parsing blast hits
    blast_hits = load_blast_results(args.blast_hits, perc_len=args.perc_len)
    primer_info = hits2info(blast_hits,
                            min_amplicon_len=args.min_amplicon_len,
                            max_amplicon_len=args.max_amplicon_len)
    nonspec_primers = [k for k,v in primer_info.items() if v[0] == 0]
    # filtering primer fasta
    ## degenerate primers
    primers_degen = load_primer_fasta(args.primer_degen_fasta, nonspec_primers, seqtype='degen')
    write_primer_fasta(primers_degen, args.prefix, seqtype='degen')
    ## expanded primers
    primers_exp = load_primer_fasta(args.primer_expand_fasta, nonspec_primers, seqtype='expand')
    write_primer_fasta(primers_exp, args.prefix, seqtype='expand')
    ## primer metadata
    edit_primer_info(args.primer_info, primer_info, primers_degen, args.prefix)
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
