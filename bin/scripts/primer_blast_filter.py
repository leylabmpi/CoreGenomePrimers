#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import argparse
import logging
import collections

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

# argparse
desc = 'Filtering primer blast output'
epi = """DESCRIPTION:
Filtering primer blast output and compiling information about
passed primers
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('blast_results', metavar='blast_results', type=str,
                    help='Blast results (-outfmt=6)')
parser.add_argument('primer_degen_fasta', metavar='primer_degen_fasta', type=str,
                    help='Degenerate primer fasta')
parser.add_argument('primer_expand_fasta', metavar='primer_expand_fasta', type=str,
                    help='Expanded-degens primer fasta')
parser.add_argument('--primer-info', type=str, default=None,
                    help='Primer info file generated by primerprospector')
parser.add_argument('--target-fasta', type=str, default=None,
                    help='Fasta of target sequences. blast hits to these will be ignored')
parser.add_argument('--perc-len', type=float, default=80,
                    help='Hit is not legit if (aln_len/query_len*100) < X')
parser.add_argument('--prefix', type=str, default='blast_filtered_primers',
                    help='Output file prefix')
parser.add_argument('--min-amplicon-len', type=int, default=100,
                    help='Min. amplicon size to be considered an off-target hit')
parser.add_argument('--max-amplicon-len', type=int, default=1500,
                    help='Max. amplicon size to be considered an off-target hit')
parser.add_argument('--version', action='version', version='0.0.1')

# functions
def load_target_fasta(infile):
    """
    Reading in target fasta file
    """
    if infile is None:
        return set()
    logging.info('Reading file: {}'.format(infile))
    gene_ids = set()
    with open(infile) as inF:
        for i,line in enumerate(inF):
            line = line.rstrip()
            if line == '':
                continue
            if line.startswith('>'):
                gene_ids.add(line.lstrip('>'))
    logging.info('  No. of target gene seqs: {}'.format(len(gene_ids)))
    return gene_ids

def to_pos_strand(sstart, send, slen):
    """
    If blast position info is scaled on the negative strand, flipping the info to the positive
    """
    # nothing required if positive strand
    if send >= sstart:
        return sstart,send,'+'
    # convert to positive strand position if on negative
    return slen-sstart,slen-send,'-'

def load_blast_results(infile, perc_len=90.0, targets=set()):
    """
    Expected input format: [qaccver saccver pident length mismatch gapopen qstart 
                            qend sstart send evalue bitscore slen qlen]
    
    Return: {primer_id : {fwd|rev : [[saccver, sstart, send],...]}}
    """
    logging.info('Reading file: {}'.format(infile))
    regex = re.compile(r'^([0-9]+)([fri])_([0-9]+)$')
    status = {'non-target' : 0, 'target' : 0, 'non-sig' : 0}
    blast_hits = collections.defaultdict(dict)
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            # filter by targets (OK if blast hit)
            if line[1] in targets:
                status['target'] += 1
                continue
            # filter by perc. alignment length
            aln_len = float(line[3]) / float(line[-1]) * 100
            if aln_len < perc_len:
                status['non-sig'] += 1
                continue
            # keeping primer
            try:
                primer_id,fri,derep_id = regex.search(line[0]).groups()
            except AttributeError:
                msg = 'Primer name formatted incorrectly: {}'
                raise AttributeError(msg.format(line[0]))
            saccver = line[1]
            slen = int(line[-2])             
            sstart = int(line[8])
            send = int(line[9])
            sstart,send,strand = to_pos_strand(sstart, send, slen)
            try:
                blast_hits[primer_id][saccver][fri].append([sstart, send, strand])
            except KeyError:
                try:
                    blast_hits[primer_id][saccver][fri] = [[sstart, send, strand]]
                except KeyError:
                    blast_hits[primer_id][saccver] = {fri : [[sstart, send, strand]]}
            status['non-target'] += 1
    # status
    if len(targets) > 0:
        logging.info('  No. of target hits: {}'.format(status['target']))
        logging.info('  No. of non-target hits: {}'.format(status['non-target']))
    else:
        cnt = status['non-target'] + status['target']
        logging.info('  No. of signfiicant hits: {}'.format(cnt))
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

def filter_blast_by_amplicon(blast_hits, min_amplicon_len, max_amplicon_len):
    """
    Filtering primers by putative amplicon that would be generated.
    If the amplicon size is outsize of the min/max, then the primers not legit off-targets.
    """
    logging.info('Filtering to only hits producing a legitimate amplicon...')
    nonspec_primers = set()
    for primer_id,d in blast_hits.items():
        status = {'no_amp' : 0, 'hit' : 0, 'wrong_strand' : 0}
        for saccver,dd in d.items():
            if primer_id in nonspec_primers:
                break
            # hits for primer pair?
            try:
                _,_ = dd['f'],dd['r']
            except KeyError:
                status['no_amp'] += 1
                continue
            # calc amplicon size of any expanded fwd-rev pair
            for x in dd['f']:
                if primer_id in nonspec_primers:
                    break
                for y in dd['r']:
                    amp_len = calc_amp_len(x[0], x[1], y[0], y[1])
                    if (x[2] != y[2] and amp_len >= min_amplicon_len
                        and amp_len <= max_amplicon_len):
                        # legit hit: different strand & amplicon_len w/in size range
                        nonspec_primers.add(primer_id)
                        status['hit'] += 1
                        break
                    elif (x[2] == y[2] and amp_len >= min_amplicon_len
                          and amp_len <= max_amplicon_len):
                        # same strand, but correct amplicon size
                        status['wrong_strand'] += 1
        # summary
        msg = '  Primer {}: legit amplicon: {}, no amplicon: {}'
        logging.info(msg.format(primer_id, status['hit'], status['no_amp']))
    # summary
    msg = '  No. of primers producing a legit non-target amplicon: {}'
    logging.info(msg.format(len(nonspec_primers)))
    return nonspec_primers

def load_primer_fasta(infile, nonspec_primers, seqtype='degen'):
    """
    Loading primer fasta. Excluding non-specific primers (nonspec_primers).    
    """
    logging.info('Reading file: {}'.format(infile))
    seqs = collections.defaultdict(dict)
    regex = re.compile(r'^[0-9]+([fri])')
    seqid = None
    n_seq = 0
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            if line.startswith('>'):
                seqid = line.lstrip('>')
                try:
                    fri = regex.search(seqid).groups()[0]
                except AttributeError:
                    msg = 'Cannot find "^[0-9]+([fri])" pattern in "{}"'
                    raise AttributeError(msg.format(seqid))
                n_seq += 1
                continue
            # primer id 
            primer_id = seqid.split('_')[0].rstrip(fri)
            if seqtype == 'degen':
                num = 0
            elif seqtype == 'expand':
                num = int(seqid.split('_')[1])
            else:
                raise ValueError('seqtype not recognized')
            # adding sequences
            try:
                seqs[primer_id][fri][num] += line
            except KeyError:
                try:
                    seqs[primer_id][fri][num] = line
                except KeyError:
                    seqs[primer_id][fri] = {num : line}
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
    outfile = prefix + '_{}.fna'.format(seqtype)
    cnt = 0
    with open(outfile, 'w') as outF:
        for primer,d in primers.items():
            for fri,dd in d.items():
                for num in sorted(dd.keys()):
                    if seqtype == 'degen':
                        outF.write('>{}{}\n{}\n'.format(primer, fri, d[fri][num]))
                    elif seqtype == 'expand':
                        outF.write('>{}{}_{}\n{}\n'.format(primer, fri, num, d[fri][num]))
                    else:
                        raise ValueError('seqtype not recognized')
                    cnt += 1
    logging.info('File written: {}'.format(outfile))
    logging.info('  No. of primer pairs written: {}'.format(cnt))

def edit_primer_info(infile, primers, prefix):
    """
    Loading info on primers, editing, and writing
    """
    if infile is None:
        return None
    outfile = prefix + '.tsv'
    
    logging.info('Reading file: {}'.format(infile))
    header = {}
    cnt = 0
    with open(infile) as inF, open(outfile, 'w') as outF:
        for i,line in enumerate(inF):
            line = line.rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                outF.write('\t'.join(line) + '\n')
                continue
            primer_id = str(line[header['primer_set']])
            try:
                x = primers[primer_id]
            except KeyError:
                continue
            if len(x.keys()) > 0:
                outF.write('\t'.join(line) + '\n')
                cnt += 1
    logging.info('  File written: {}'.format(outfile))            
    logging.info('  No. of records written: {}'.format(cnt))

def main(args):
    """
    Main interface
    """
    # output
    outdir = os.path.split(args.prefix)[0]
    if outdir != '' and not os.path.isdir(outdir):
        os.makedirs(outdir)    
    # target gene IDs
    target_gene_ids = load_target_fasta(args.target_fasta)    
    # which primers hit non-targets (true target hits ignored)?
    blast_hits = load_blast_results(args.blast_results, args.perc_len, targets=target_gene_ids)
    nonspec_primers = filter_blast_by_amplicon(blast_hits,
                                               args.min_amplicon_len,
                                               args.max_amplicon_len)
    # filtering primers
    ## degenerate primers
    primers_degen = load_primer_fasta(args.primer_degen_fasta, nonspec_primers, seqtype='degen')
    write_primer_fasta(primers_degen, args.prefix, seqtype='degen')
    ## expanded primers
    primers_exp = load_primer_fasta(args.primer_expand_fasta, nonspec_primers, seqtype='expand')
    write_primer_fasta(primers_exp, args.prefix, seqtype='expand')
    ## primer metadata
    edit_primer_info(args.primer_info, primers_degen, args.prefix)
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)