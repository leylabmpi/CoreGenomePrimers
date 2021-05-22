#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import gzip
import uuid
import argparse
import logging
import collections

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Filtering two fasta files down to just intersection, formatting, & writting metadata'
epi = """DESCRIPTION:
Filtering 2 fasta files down to the intersection of their sequencines.
The sequence headers must perfectly match.
If any duplicate headers, only the first will be selected.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('nuc_fasta', metavar='nuc_fasta', type=str,
                    help='nucleotide fasta file')
parser.add_argument('prot_fasta', metavar='prot_fasta', type=str,
                    help='amino acid fasta file')
parser.add_argument('gff', metavar='gff', type=str,
                    help='gff file')
parser.add_argument('taxon', metavar='taxon', type=str,
                    help='taxon name')
parser.add_argument('--nuc-fasta-out', type=str, default='formatted.fna',
                    help='Reformatted nucleotide cds fasta')
parser.add_argument('--prot-fasta-out', type=str, default='formatted.faa',
                    help='Reformatted amino acid cds fasta')
parser.add_argument('--metadata-out', type=str, default='formatted.tsv',
                    help='gene metadata (tab-delim format)')
parser.add_argument('--rrna-fasta-out', type=str, default='formatted.ffn',
                    help='Reformatted nucleotide rRNA fasta')
parser.add_argument('--version', action='version', version='0.0.1')

def make_index(fasta):
    """
    Creating set of fasta names.
    Parsing by CDS and rRNA sequence IDs
    """
    regex1 = re.compile(r' .+')
    idx = set()
    with open(fasta) as inF:
        for line in inF:
            if line.startswith('>'):
                line = line.lstrip('>').rstrip()
                idx.add(regex1.sub('', line))
    print('KOPKOCEE_00422' in idx)
    return idx

def filter_fasta(fasta, idx, output, gzip_out=False, rRNA=False):
    """
    Filtering fasta to just those in idx
    """
    found = {}
    hit = False
    regex= re.compile(r' .+')
    regex_rrna = re.compile(r' [0-9]+S ribosomal RNA')
    with open(fasta) as inF, open(output, 'w') as outF:
        for line in inF:
            if fasta.endswith('.gz'):
                line = line.decode('utf8')
            if line.startswith('>'):
                is_rRNA = regex_rrna.search(line) is not None
                line = regex.sub('', line.lstrip('>').rstrip())
                # filter is already seen
                try:
                    found[line]
                    continue
                except KeyError:
                    pass
                # is rRNA?
                if is_rRNA is True:
                    if rRNA is False:
                        hit = False
                        continue
                else:
                    if rRNA is True:
                        hit = False
                        continue
                # is seq in index?
                try:                    
                    found[line] = idx[line]
                    hit = True
                except KeyError:
                    hit = False
                    continue
                # writing
                seq_name = '>' + idx[line] + '\n'
                outF.write(seq_name)
            else:
                if hit:
                    outF.write(line)
                        
    logging.info('File written: {}'.format(output))
    logging.info('Number of seqs written: {}'.format(len(found.keys())))

def write_meta(infile, idx, taxon, outfile):
    """
    Writing metadata file
    """
    logging.info('Reading file: {}'.format(infile))
    header = ['seq_uuid', 'seq_orig_name', 'contig_id',
              'taxon',  'start', 'end',
              'score', 'strand', 'annotation']
    seq_cnt = set()
    with open(infile) as inF, open(outfile, 'w') as outF:
        outF.write('\t'.join(header) + '\n')
        for i,line in enumerate(inF):
            if line.startswith('##FASTA'):
                break             
            if line.startswith('#'):
                continue
            line = line.rstrip().split('\t')
            # seqid
            line_idx = {}
            for x in line[-1].split(';'):
                y = x.split('=')
                line_idx[y[0]] = y[1]
            try:
                seq_id = line_idx['ID']
            except KeyError:
                msg = 'Cannot find "ID" tag in Line {}'
                logging.warning(msg.format(i+1))
                continue
            try:
                annotation = line_idx['product']
            except KeyError:
                msg = 'Cannot find "product" tag in Line {}'
                logging.warning(msg.format(i+1))
                annotation = 'not annotated'
            # getting uuid
            try:
                seq_uuid = idx[seq_id]
            except KeyError:
                msg = 'Cannot find sequence id in the index: {}'
                raise KeyError(msg.format(seq_id))
            # output line
            line = [str(x) for x in line]
            line = [seq_uuid, seq_id, line[0], taxon] + line[3:-2] + [annotation]
            outF.write('\t'.join(line) + '\n')
            seq_cnt.add(seq_uuid)
    # status
    logging.info('File written: {}'.format(outfile))
    logging.info('Number of records written: {}'.format(len(seq_cnt)))
    
def main(args):
    """
    Main interface
    """    
    # creating the seq header index
    seq_idx = make_index(args.nuc_fasta)
    seq_idx = {x:str(uuid.uuid4()).replace('-', '') for x in seq_idx}
    # formatting the fasta files
    filter_fasta(args.nuc_fasta, seq_idx, args.nuc_fasta_out)
    filter_fasta(args.prot_fasta, seq_idx, args.prot_fasta_out)
    filter_fasta(args.nuc_fasta, seq_idx, args.rrna_fasta_out, rRNA=True)     
    # creating name index
    write_meta(args.gff, seq_idx, args.taxon, args.metadata_out)
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
