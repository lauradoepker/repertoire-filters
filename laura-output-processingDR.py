#!/usr/bin/env python
import csv
import sys
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse
import numpy

partis_path = '.'  # edit this if you're not running from the main partis dir
sys.path.insert(1, partis_path + '/python')
import utils
import glutils
from clusterpath import ClusterPath

parser = argparse.ArgumentParser()
parser.add_argument('--infile')
parser.add_argument('--locus')
parser.add_argument('--param')
args = parser.parse_args()

glfo = glutils.read_glfo(args.param + '/hmm/germline-sets', locus=args.locus)

cp = ClusterPath()
cp.readfile(args.infile)
best_partition = cp.partitions[cp.i_best]
sorted_clusters = sorted(best_partition, key=len, reverse=True)  # sort by size

biggest_cluster = sorted_clusters[0]  # let's look for the biggest one in the annotation file, just for kicks
ten_biggest = sorted_clusters[0:10]  # let's look at the 10 biggest clusters
#for cluster in ten_biggest:
    #cluster_index = ten_biggest.index(cluster)

bnab_genes = ['IGHV1-2', 'IGHV1-46', 'IGHV1-3', 'IGHV1-69', 'IGHV3-23', 'IGHV3-30', 'IGHV4-59']

def print_stuff(line):
    print line['v_gene'],
    print line['d_gene'],
    print line['j_gene']
    print len(line['unique_ids']), 'sequences in cluster'
    print 'number of mutations per sequence in cluster', sorted(line['n_mutations'])
    print numpy.mean(line['n_mutations']), 'mean number of mutations'
    print numpy.median(line['n_mutations']), 'median number of mutations'
    print len(line['naive_seq']), 'length of naive seq'
    print numpy.mean(line['n_mutations'])/len(line['naive_seq']), 'mean % SHM'
    print line['cdr3_length'], 'CDR3 length in nts'
    print

def getkey(uid_list):
    return ':'.join(uid_list)

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
def keyfunc(q):
    return len(annotations[q]['unique_ids'])

sorted_clusters = sorted(annotations, key=keyfunc, reverse=True)

# # loop over ten biggest
# for cluster in sorted_clusters[:10]:
#     print len(annotations[cluster]['unique_ids'])

# add more criteria
def boolfunc(q):
    if annotations[cluster]['cdr3_length'] < 50:
        return False
    if len(annotations[cluster]['unique_ids']) < 20:
        return False
    return True

interesting_clusters = [cluster for cluster in sorted_clusters if boolfunc(cluster)]
print '  found %d interesting clusters' % len(interesting_clusters)
for cluster in interesting_clusters:
    print_stuff(annotations[cluster])

sys.exit()

if line['unique_ids'] in ten_biggest:
    print 'found the 10 biggest clusters!'
    #print line['n_mutations']
    print line['v_gene'],
    print line['d_gene'],
    print line['j_gene']
    print len(line['unique_ids']), 'sequences in cluster'
    print 'number of mutations per sequence in cluster', sorted(line['n_mutations'])
    print numpy.mean(line['n_mutations']), 'mean number of mutations'
    print numpy.median(line['n_mutations']), 'median number of mutations'
    print len(line['naive_seq']), 'length of naive seq'
    print numpy.mean(line['n_mutations'])/len(line['naive_seq']), 'mean % SHM'
    print line['cdr3_length'], 'CDR3 length in nts'
    print
    # for key, val in line.items():
        # print '%20s %s' % (key, val)

    # utils.print_reco_event(line) # print ascii-art representation of the rearrangement event

if line['cdr3_length'] > 80 and len(line['unique_ids']) > 20: #large family with long cdr3
    print 'found big clusters with long cdr3 lengths!'
    print line['v_gene'],
    print line['d_gene'],
    print line['j_gene']
    print len(line['unique_ids']), 'sequences in cluster'
    print 'number of mutations per sequence in cluster', sorted(line['n_mutations'])
    print numpy.mean(line['n_mutations']), 'mean number of mutations'
    print numpy.median(line['n_mutations']), 'median number of mutations'
    print len(line['naive_seq']), 'length of naive seq'
    print numpy.mean(line['n_mutations'])/len(line['naive_seq']), 'mean % SHM'
    print line['cdr3_length'], 'CDR3 length in nts'
    print
    # for key, val in line.items():
        # print '%20s %s' % (key, val)
    # utils.print_reco_event(line) # print ascii-art representation of the rearrangement event

if (line['v_gene']).split('*')[0] in bnab_genes and len(line['unique_ids']) > 50: #large family with known bnAb VH gene
    print 'found big clusters with bnAb gene usage!'
    print line['v_gene'],
    print line['d_gene'],
    print line['j_gene']
    print len(line['unique_ids']), 'sequences in cluster'
    print 'number of mutations per sequence in cluster', sorted(line['n_mutations'])
    print numpy.mean(line['n_mutations']), 'mean number of mutations'
    print numpy.median(line['n_mutations']), 'median number of mutations'
    print len(line['naive_seq']), 'length of naive seq'
    print numpy.mean(line['n_mutations'])/len(line['naive_seq']), 'mean % SHM'
    print line['cdr3_length'], 'CDR3 length in nts'
    print
    # for key, val in line.items():
        # print '%20s %s' % (key, val)
    # utils.print_reco_event(line) # print ascii-art representation of the rearrangement event
