#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import pickle
import uuid
import collections
import argparse
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Formatting blastx results'
epi = """DESCRIPTION:
* input format: qaccver saccver pident length mismatch qstart qend sstart send evalue slen qlen sscinames staxids
* formatting raw blastx results
  * eg., including gene names for each accession
* tsv file written (with header)
  * output gzip'ed
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('blast_hits', metavar='blast_hits', type=str,
                    help='BLAST hits file')
parser.add_argument('index', metavar='index', type=str,
                    help='NCBI-nr Accession<->GeneName index')
parser.add_argument('outfile', metavar='outfile', type=str,
                    help='Output file name')
parser.add_argument('--taxids', type=str, default=None,
                    help='File of taxids used for filtering the blast table')
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

def load_index(infile):
    logging.info('Reading file: {}'.format(infile))
    idx = {}    
    if os.path.splitext(infile)[1] == '.pkl':
        # if pickled file
        logging.info('  Assuming the file is pickled...')
        with open(infile, 'rb') as inF:
            idx = pickle.load(inF)
    else:
        # else assuming tsv file
        with _open(infile) as inF:
            for line in inF:
                line = _decode(line).rstrip().split('\t')
                if line[0] == 0:
                    continue
                idx[line[0]] = line[1]
    # status
    logging.info('  No. of accessions: {}'.format(len(idx.keys())))
    return idx

def load_taxids(infile):
    """
    Loading taxids; assuming 1 per line
    """
    if infile is None:
        return None
    taxids = set()
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            taxids.add(int(line[0]))
    return taxids

def filter_by_taxid(line, taxids):
    if taxids is None:
        return False
    for x in line[-1].split(';'):
        if int(x) in taxids:
            msg = 'taxid "{}" in --taxids, filtering record: [{},{}]'
            logging.warning(msg.format(x, line[0], line[1]))
            return True
    return False

def main(args):
    # loading taxids
    taxids = load_taxids(args.taxids)
    
    # loading index
    idx = load_index(args.index)

    # reading in blast results & formatting
    logging.info('Reading file: {}'.format(args.blast_hits))
    hits = {}    
    with open(args.blast_hits) as inF:
        for line in inF:
            line = line.rstrip().split('\t')
            # filtering by taxid
            if filter_by_taxid(line, taxids) is True:
                continue
            # getting gene annotation
            seqid,cluster_id = line[0].split('_', 1)
            try:
                gene_name = idx[line[1]]
            except KeyError:
                logging.warning('Cannot find accession in index: {}'.format(line[1]))
                gene_name = 'None'
            # saving output
            try:
                hits[cluster_id].append([seqid, line[1], gene_name] + line[2:])
            except KeyError:
                hits[cluster_id] = [[seqid, line[1], gene_name] + line[2:]]
                
    # ranking hits by pident
    for cluster_id,records in hits.items():
        pident = [x[3] for x in records]
        ranks = [sorted(pident, reverse=True).index(x)+1 for x in pident]
        for record,rank in zip(records,ranks):
            record.append(str(rank))

    # writing hits
    regex = re.compile(r'^cluster_')
    header = ['cluster_id', 'query', 'subject', 'subject_name',
              'pident', 'length', 'mismatch', 'qstart', 'qend',
              'sstart', 'send', 'evalue', 'slen', 'qlen', 'sscinames',
              'staxids', 'pident_rank']    
    with open(args.outfile, 'w') as outF:
        # header
        outF.write('\t'.join(header) + '\n')
        # body
        for cluster_id,records in hits.items():            
            for record in records:
                outF.write('\t'.join([regex.sub('', cluster_id)] + record) + '\n')
    logging.info('File written: {}'.format(args.outfile))

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
