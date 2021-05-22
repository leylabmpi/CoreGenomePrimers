#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import argparse
import logging
import collections
from pyfaidx import Fasta

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Determining which gene clusters are core genes'
epi = """DESCRIPTION:
A list of IDs for core clusters (one per line) written to STDOUT
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('membership', metavar='membership', type=str,
                    help='mmseqs cluster membership file')
parser.add_argument('uclust', metavar='uclust', type=str,
                    help='uclust file')
parser.add_argument('metadata', metavar='metadata', type=str,
                    help='genes metadata file')
parser.add_argument('nuc_fasta', metavar='nuc_fasta', type=str,
                    help='nucleotide fasta file')
parser.add_argument('aa_fasta', metavar='aa_fasta', type=str,
                    help='nucleotide fasta file')
parser.add_argument('--outdir', type=str, default='clusters',
                    help='Output directory')
parser.add_argument('--perc-genomes-rrna', type=float, default=100.0,
                    help='Gene must be found in >=X percent of genomes')
parser.add_argument('--perc-genomes-cds', type=float, default=100.0,
                    help='Gene must be found in >=X percent of genomes')
parser.add_argument('--copies-per-genome-rrna', type=int, default=10,
                    help='Max gene copy per genome')
parser.add_argument('--copies-per-genome-cds', type=int, default=1,
                    help='Max gene copy per genome')
parser.add_argument('--max-clusters-rrna', type=int, default=0,
                    help='Max number of rrna clusters to write out. If <1, all written.')
parser.add_argument('--max-clusters-cds', type=int, default=0,
                    help='Max number of cds clusters to write out. If <1, all written.')
parser.add_argument('--version', action='version', version='0.0.1')

def load_meta(infile):
    """
    return: {seq_uuid : taxon_id}
    """
    logging.info('Loading file: {}'.format(infile))
    header = {}
    idx = {}
    taxon_idx = {}
    with open(infile) as inF:
        for i,line in enumerate(inF):
            line = line.rstrip().split('\t')        
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                continue
            if line[0] == '':
                continue
            # converting taxon to number
            taxon = line[header['taxon']]
            try:
                cnt = taxon_idx[taxon]
            except KeyError:
                cnt = len(taxon_idx.keys()) + 1
                taxon_idx[taxon] = cnt        
            idx[line[header['seq_uuid']]] = cnt
    logging.info('No. of seqs loaded: {}'.format(len(idx.keys())))
    return idx, len(set(idx.values()))

def load_uclust(infile, idx):
    """
    return: {cluster : {taxon : [gene_count, 'rrna|cds'}}
    """
    logging.info('Loading file: {}'.format(infile))
    clst = collections.defaultdict(dict)
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            # getting taxon for seq in cluster
            try:
                taxon = idx[line[8]]
            except KeyError:
                msg = 'Cannot find taxon for {}'
                raise KeyError(msg.format(line[8]))
            # counting
            try:
                clst[line[1]][taxon] += 1                
            except KeyError:
                clst[line[1]][taxon] = 1
    logging.info('  No. of clusters loaded: {}'.format(len(clst.keys())))
    return clst

def load_membership(infile, idx):
    """
    return: {cluster : {taxon : gene_count}}
    """
    logging.info('Loading file: {}'.format(infile))
    clst = collections.defaultdict(dict)
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            # getting taxon for seq in cluster
            try:
                taxon = idx[line[1]]
            except KeyError:
                msg = 'Cannot find taxon for {}'
                raise KeyError(msg.format(line[1]))
            # counting
            try:
                clst[line[0]][taxon] += 1
            except KeyError:
                clst[line[0]][taxon] = 1
    idx.clear()
    logging.info('  No. of clusters: {}'.format(len(clst.keys())))
    return clst

def is_core(d, ntaxa, frac, copies_per_genome):
    """
    True if present for all taxa & one copy per each taxon
    """
    # fraction of taxa that gene is present
    if round(len(d.keys()) / float(ntaxa),2) < frac:
        return False
    # copy number
    if any([x > copies_per_genome for x in d.values()]):
        return False
    return True

def get_core_clusters(clst, genetype, ntaxa, perc_genomes, max_clusters, copies_per_genome):
    """
    Which gene clusters are core?
    """
    frac = perc_genomes / 100.0
    logging.info('Finding core {} clusters...'.format(genetype))
    core_clst = set()
    for cluster,d in clst.items():
        if is_core(d, ntaxa, frac, copies_per_genome) is True:
            core_clst.add(cluster)         
    clst.clear()
    msg = '  No. of core clusters: {}'
    logging.info(msg.format(len(core_clst)))
    if max_clusters >= 1:
        core_clst = set(list(core_clst)[:max_clusters])
        logging.info('  Only keeping {} clusters'.format(max_clusters))
    return core_clst                   

def get_seq_ids(infile, core):
    """
    Getting seqIDs for each gene cluster
    """
    logging.info('Getting seq IDs for core clusters...')
    logging.info('  Loading file: {}'.format(infile))
    idx = {}
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            # {clusterid : seqid}
            if line[0] in core:
                try:
                    idx[line[0]].append(line[1])
                except KeyError:
                    idx[line[0]] = [line[1]]
    logging.info('  No. of clusters: {}'.format(len(idx.keys())))
    return idx

def get_seq_ids_uclust(infile, core):
    """
    Getting seqIDs for each gene cluster
    """
    logging.info('Getting seq IDs for core clusters...')
    logging.info('  Loading file: {}'.format(infile))
    idx = {}
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            # {clusterid : seqid}
            if line[1] in core:
                try:
                    idx[line[1]].append(line[8])
                except KeyError:
                    idx[line[1]] = [line[8]]
    logging.info('  No. of clusters: {}'.format(len(idx.keys())))
    return idx

def parse_seqs(genes, core, outdir, ext, allow_missing=False):
    """
    Writing out gene clusters
    """
    logging.info('Writing out {} sequences...'.format(ext))
    clust_idx = collections.defaultdict(dict)
    for genetype in core.keys():
        D = os.path.join(outdir, genetype)
        if not os.path.isdir(D):
            os.makedirs(D)
        n_files = 0
        for i,(cluster,seqs) in enumerate(core[genetype].items()):
            i += 1
            clust_idx[genetype][cluster] = i
            outfile = os.path.join(D, 'cluster_{}.{}'.format(i, ext))
            with open(outfile, 'w') as outF:
                for seq in seqs:
                    try:
                        outF.write('>' + genes[seq].name + '\n')
                        outF.write(str(genes[seq]) + '\n')
                    except KeyError:
                        if allow_missing:
                            continue
                        else:
                            msg = 'Cannot find {} in gene sequences'
                            raise KeyError(msg.format(seq))
                n_files += 1
            if os.stat(outfile).st_size == 0:
                n_files -= 1
                os.unlink(outfile)
        logging.info('  No. of {} fasta files written: {}'.format(genetype, n_files))
    return clust_idx

def clust2metadata(core, metadata_file, membership_file, clust_idx, outdir):
    logging.info('Formatting core gene metadata table...')
    # creating {genetype : {core_seq_id : core_cluster_id}}
    idx = collections.defaultdict(dict)
    for genetype in core.keys():
        for cluster,seqs in core[genetype].items():
            for seq in seqs:
                idx[genetype][seq] = cluster
    core.clear()
    n_core = sum([len(v.values()) for v in idx.values()])
    logging.info('  No. of core seqs: {}'.format(n_core))
    # loading gene records by cluster
    clusts = collections.defaultdict(dict)
    header = collections.OrderedDict()
    with open(metadata_file) as inF:
        for i,line in enumerate(inF):
            line = line.rstrip().split('\t')
            if i == 0:
                for ii,x in enumerate(line):
                    header[x] = ii
                continue
            for genetype in idx.keys():
                try:
                    clust_name = idx[genetype][line[header['seq_uuid']]]                
                except KeyError:
                    continue
                try:
                    clusts[genetype][clust_name].append(line)
                except KeyError:
                    clusts[genetype][clust_name] = [line]
    ## status
    for genetype in clusts.keys():
        n_clusts = len(clusts[genetype].keys())
        logging.info('  No. of {} clusters: {}'.format(genetype, n_clusts))
        nseqs = sum([len(x) for x in clusts[genetype].values()])
        logging.info('  No. of {} gene records: {}'.format(genetype, nseqs))
    # writing files
    for genetype in clusts.keys():
        D = os.path.join(outdir, genetype)
        if not os.path.isdir(D):
            os.makedirs(D)
        cnt = 0
        for clust_name,lines in clusts[genetype].items():
            clust_id = clust_idx[genetype][clust_name]
            outfile = os.path.join(D, 'cluster_{}.tsv'.format(clust_id))
            with open(outfile, 'w') as outF:
                cnt += 1
                X = [str(x) for x in header.keys()] + ['cluster_name', 'clust_id', 'gene_type']
                outF.write('\t'.join(X) + '\n')
                for line in lines:
                    Y = [str(clust_name), str(clust_id), genetype]
                    outF.write('\t'.join(line + Y) + '\n')
        logging.info('  No. of {} gene metadata files written: {}'.format(genetype, cnt))
    
def main(args):
    """
    Main interface
    """
    if args.outdir != '' and not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    # load metdata {seq_uuid : taxon}
    idx, ntaxa = load_meta(args.metadata)
    # load clusters
    clst = {}
    clst['rrna'] = load_uclust(args.uclust, idx)
    clst['cds'] = load_membership(args.membership, idx)
    # finding core clusters
    core = {}
    core['rrna'] = get_core_clusters(clst['rrna'], 'rrna', ntaxa,
                                     args.perc_genomes_rrna,
                                     args.max_clusters_rrna,
                                     args.copies_per_genome_rrna)
    core['cds'] = get_core_clusters(clst['cds'], 'cds', ntaxa,
                                    args.perc_genomes_cds,
                                    args.max_clusters_cds,
                                    args.copies_per_genome_cds)
    # getting seq ids per core cluster
    core_ids = {}
    core_ids['rrna'] = get_seq_ids_uclust(args.uclust, core['rrna'])
    core_ids['cds'] = get_seq_ids(args.membership, core['cds'])
    # parsing sequence    
    ## indexing, then writing sequence files per cluster
    genes = Fasta(args.nuc_fasta)
    clust_idx = parse_seqs(genes, core_ids, args.outdir, ext='fna')
    ## indexing, then writing sequence files per cluster
    genes = Fasta(args.aa_fasta)
    clust_idx = parse_seqs(genes, core_ids, args.outdir, ext='faa', allow_missing=True)
    # core gene cluster metadata with cluster info
    clust2metadata(core_ids, args.metadata, args.membership, clust_idx, args.outdir)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
