#!/usr/bin/env python
import csv
import sys
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse
import numpy
from Bio.Seq import Seq

partis_path = '.'  # edit this if you're not running from the main partis dir
sys.path.insert(1, partis_path + '/python')
import utils
import glutils
from clusterpath import ClusterPath

parser = argparse.ArgumentParser()
parser.add_argument('--infile')
parser.add_argument('--locus')
parser.add_argument('--param')
parser.add_argument('--cdr3')
args = parser.parse_args()

glfo = glutils.read_glfo(args.param + '/hmm/germline-sets', locus=args.locus)
print(sys.argv)
print 'infile =', args.infile
print 'param =', args.param
print

cp = ClusterPath()
cp.readfile(args.infile)
best_partition = cp.partitions[cp.i_best]

# clonal family attributes to print
def print_stuff(line):
    cluster_index = sorted_clusters.index(cluster)
    naive_cdr3, matureiseq0_cdr3 = utils.subset_sequences(line, iseq=0, restrict_to_region='cdr3') # returns the CDR3 nt sequence for naive, and the first mutated sequence (iseq0); CDR3 = first base of cysteine through last base of tryptophan

    # mature_cdr3_seqs = []  # trying to translate the consensus cdr3 so I can search these with my seed seqs
    # for iseq in range(len(line['unique_ids'])):
    #     naive_cdr3_seq, mature_cdr3_seq = utils.subset_sequences(line, iseq=iseq, restrict_to_region='cdr3')
    #     mature_cdr3_seqs.append(mature_cdr3_seq)
    # mature_cdr3_seqs
    # translated_cdr3 = mature_cdr3_seqs.translate()

    cdr3_aa = '%-30s' % Seq(naive_cdr3).translate()
    # If a cluster contains one of our seed seqs, color this CDR3 red
    if any('-ig' in s for s in line['unique_ids']):
        cdr3_aa = utils.color('red', cdr3_aa, width=30)
    if args.cdr3 in cdr3_aa: # Only print clusters with naive CDR3 that matches our specified --cdr3 argument
        print 'index    genes                                        size    n muts    SHM     rep frac     CDR3                                FayWuH'
        print '                                                            mean  med                        len  seq'
        print '%4s     %s %s %s %5d %5d %5d %7.3f   %8.4f     %2d   %s %4.2f' % (
                cluster_index,
                utils.color_gene(line['v_gene'], width=15),
                utils.color_gene(line['d_gene'], width=15),
                utils.color_gene(line['j_gene'], width=10),
                len(line['unique_ids']),
                numpy.mean(line['n_mutations']),
                numpy.median(line['n_mutations']),
                numpy.mean(line['mut_freqs']),
                float(len(cluster)) / n_total,
                (line['cdr3_length']/3),
                cdr3_aa,
                utils.fay_wu_h(line, debug=False),
                )
        # print 'number of mutations per sequence in cluster', sorted(line['n_mutations'])
        print len(line['naive_seq']), 'length of naive seq'
        # utils.print_reco_event(utils.synthesize_single_seq_line(line, iseq=0))  # print ascii-art representation of the rearrangement event
        print 'unique_ids: ', getkey(line['unique_ids'])
        print
        print utils.print_reco_event(line)

# formatting necessity
def getkey(uid_list):
    return ':'.join(uid_list)

# creates a dictionary with keys = unique_ids and values = annotations
annotations = {}
with open(args.infile.replace('.csv', '-cluster-annotations.csv')) as csvfile:
    reader = csv.DictReader(csvfile)
    for line in reader:  # there's a line for each cluster
        if line['v_gene'] == '':  # failed (i.e. couldn't find an annotation)
            continue
        utils.process_input_line(line)  # converts strings in the csv file to floats/ints/dicts/etc.
        utils.add_implicit_info(glfo, line)  # add stuff to <line> that's useful, isn't written to the csv since it's redundant
        # utils.print_reco_event(line)  # print ascii-art representation of the rearrangement event
        annotations[getkey(line['unique_ids'])] = line


# sort by size
sorted_clusters = sorted(annotations, key=lambda q: len(annotations[q]['unique_ids']), reverse=True)
n_total = sum([len(cluster) for cluster in sorted_clusters])

# # add more criteria
# biggest_clusters = sorted_clusters[:100] # 100 biggest clusters
# shm_clusters = sorted(biggest_clusters, key=lambda q: numpy.mean(annotations[q]['mut_freqs']), reverse=True)
# sfs_clusters = sorted(biggest_clusters, key=lambda q: utils.fay_wu_h(annotations[q], debug=False))

# cluster size: print x biggest clusters
print 'printing chosen cluster: ' + str(args.cdr3)
for cluster in sorted_clusters[:100]:
    print_stuff(annotations[cluster])

sys.exit()
