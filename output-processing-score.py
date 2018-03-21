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
parser.add_argument('--nclust')
args = parser.parse_args()

glfo = glutils.read_glfo(args.param + '/hmm/germline-sets', locus=args.locus)

print(sys.argv)
print 'infile =', args.infile
print 'param =', args.param

cp = ClusterPath()
cp.readfile(args.infile)
best_partition = cp.partitions[cp.i_best]
# sorted_clusters = sorted(best_partition, key=len, reverse=True)  # sort by size

# clonal family attributes to print
print '''

score = interest score, indicating interesting attributes: size, SHM, SFS, bnAb VH usage

Size & SHM:
4 points for rank in top 25
3 points for rank 25-50
2 points for rank 50-75
1 point for rank 75-100

SFS (Fay Wu H) scores earning 4-1 points: < -20, -10, 0, 10

4 points for using a bnAb VH gene (does not consider CDR3 length)
'''

print 'score index   genes                                       size     n muts    SHM     rep frac     CDR3                                FayWuH'
print '                                                                 mean  med                        len  seq'
def print_stuff(line):
    intscore = 0 # create a clonal family scoring system
    cluster_index = sorted_clusters.index(cluster)
    shm_index = shm_clusters.index(cluster)
    naive_cdr3, matureiseq0_cdr3 = utils.subset_sequences(line, iseq=0, restrict_to_region='cdr3') # line['naive_seq'][(line['codon_positions']['v']):((line['codon_positions']['j'])+3)] #get nt sequence of CDR3 from first base of cysteine through last base of tryptophan
    # mature_cdr3_seqs = []  # trying to translate the consensus cdr3 so I can search these with my seed seqs
    #     for iseq in range(len(line['unique_ids'])):
    #         naive_cdr3_seq, mature_cdr3_seq = utils.subset_sequences(line, iseq=iseq, restrict_to_region='cdr3')
    #         mature_cdr3_seqs.append(mature_cdr3_seq)
    # translated_cdr3 = Seq().... not done
    cdr3_aa = '%-30s' % Seq(naive_cdr3).translate()
    if any('-ig' in s for s in line['unique_ids']):
        cdr3_aa = utils.color('red', cdr3_aa, width=30)

    # score clusters based on cluster size
    if cluster_index < 25:
        intscore = intscore + 4
    elif cluster_index >= 25 and cluster_index <= 50:
        intscore = intscore + 3
    elif cluster_index >= 50 and cluster_index <= 75:
        intscore = intscore + 2
    elif cluster_index >= 75 and cluster_index <= 100:
        intscore = intscore + 1

    # score clusters based on SHM
    if shm_index < 25:
        intscore = intscore + 4
    elif shm_index >= 25 and shm_index <= 50:
        intscore = intscore + 3
    elif shm_index >= 50 and shm_index <= 75:
        intscore = intscore + 2
    elif shm_index >= 75 and shm_index <= 100:
        intscore = intscore + 1

    # score clusters based on SFS
    if utils.fay_wu_h(line, debug=False) <= -20:
        intscore = intscore + 4
    elif utils.fay_wu_h(line, debug=False) <= -10:
        intscore = intscore + 3
    elif utils.fay_wu_h(line, debug=False) <= 0:
        intscore = intscore + 2
    elif utils.fay_wu_h(line, debug=False) <= 10:
        intscore = intscore + 1

    # score by bnAb gene usage
    if (line['v_gene']).split('*')[0] in (cd4bs_genes or glycan_genes or bridging_genes or mper_genes): # beware this does not include CDR3 length of bnAb VH genes
        intscore = intscore + 4

    print '%4s %4s     %s %s %s %5d %5d %5d %7.3f   %8.4f     %2d   %s %4.2f' % (
            intscore,
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
# Troubleshooting... looking at a cluster of a certain size:
    # if len(line['unique_ids']) == 174:
    #     print line['unique_ids']
    #     print utils.print_reco_event(line)

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

#### sorted_clusters = [c for c in sorted_clusters if utils.is_functional(annotations[c])] # checks if the cluster contains ANY non-functional sequences

# total size of repertoire (number sequences)
n_total = sum([len(cluster) for cluster in sorted_clusters])

# add more criteria
biggest_clusters = sorted_clusters[:100] # 100 biggest clusters
shm_clusters = sorted(biggest_clusters, key=lambda q: numpy.mean(annotations[q]['mut_freqs']), reverse=True) # rank by SHM
sfs_clusters = sorted(biggest_clusters, key=lambda q: utils.fay_wu_h(annotations[q], debug=False)) # rank by SFS score

cluster_sfses = {}
for cluster in biggest_clusters:
    cluster_sfses[cluster] = utils.fay_wu_h(annotations[cluster], debug=False)
print numpy.mean(cluster_sfses.values())
print numpy.std(cluster_sfses.values())
print numpy.percentile(cluster_sfses.values(), 5)
print numpy.percentile(cluster_sfses.values(), 10)
print numpy.percentile(cluster_sfses.values(), 50)
print numpy.percentile(cluster_sfses.values(), 80)
print numpy.percentile(cluster_sfses.values(), 90)

# create function that gives me the score - this function calls a subfunction for each metric (i.e. percentile).  The superfunction can then weight the metrics
# give 0 points to anyone not in top 30 percentile

sys.exit()

# published bnAb VH gene usage (Yu and Guan, Frontiers in Immunology, 2014 and Wu and Kong, Curr Op in Immunol, 2016)
cd4bs_genes = ['IGHV1-2', 'IGHV1-46', 'IGHV1-3', 'IGHV4-61', 'IGHV1-69', 'IGHV3-23', 'IGHV3-30']
glycan_genes = ['IGHV3-21', 'IGHV1-8', 'IGHV3-20', 'IGHV3-33', 'IGHV4-39', 'IGHV4-59', 'IGHV4-4']
bridging_genes = ['IGHV1-3', 'IGHV3-30', 'IGHV1-28']
mper_genes = ['IGHV1-69', 'IGHV2-5', 'IGHV3-15', 'IGHV5-51']

# cluster size: print x biggest clusters
print '\x1b[1;32;40m' + '  printing the largest clusters' + '\x1b[0m'
for cluster in sorted_clusters[:5]:
    # if sorted_clusters.index(cluster) < 50:
    #     print_stuff(annotations[cluster])
    print_stuff(annotations[cluster])

# high mean %SHM: print most mutated clusters from 100 biggest clusters
mutclust = int(args.nclust)
print '\x1b[1;32;40m' + '  printing the most mutated clusters (within 100 biggest)' + '\x1b[0m'
for cluster in shm_clusters[:mutclust]:
    # if sorted_clusters.index(cluster) < 50:
    #     print_stuff(annotations[cluster])
    print_stuff(annotations[cluster])

# Highest SFS Fay Wu H scores
print '\x1b[1;32;40m' + '  printing clusters with high SFS (within 100 biggest)' + '\x1b[0m'
for cluster in sfs_clusters[:10]:
    print_stuff(annotations[cluster])


if args.locus != 'igh': #continue to bnAb VH gene usage if we're looking at heavy chain
    sys.exit()

# CD4bs bnAb VH gene usage
def cd4bs_bool(q):
    if (annotations[cluster]['v_gene']).split('*')[0] not in cd4bs_genes:
        return False
    if annotations[cluster]['cdr3_length'] > 75 or annotations[cluster]['cdr3_length'] < 30:
        return False
    return True

cd4bs_clusters = [cluster for cluster in biggest_clusters if cd4bs_bool(cluster)]
print '\x1b[1;32;40m' + '  found %d cd4bs clusters, printing the largest hits' % len(cd4bs_clusters)
sorted_cd4bs = sorted(cd4bs_clusters, key=lambda q: len(annotations[q]['unique_ids']), reverse=True)
for cluster in sorted_cd4bs[0:10]: # five biggest cd4bs clusters
    print_stuff(annotations[cluster])

# glycan epitope bnAb VH gene usage
def glycan_bool(q):
    if (annotations[cluster]['v_gene']).split('*')[0] not in glycan_genes:
        return False
    if annotations[cluster]['cdr3_length'] < 60:
        return False
    return True

glycan_clusters = [cluster for cluster in biggest_clusters if glycan_bool(cluster)]
print '\x1b[1;32;40m' + '  found %d glycan clusters, printing the largest hits' % len(glycan_clusters)
sorted_glycan = sorted(glycan_clusters, key=lambda q: len(annotations[q]['unique_ids']), reverse=True)
for cluster in sorted_glycan[0:10]: # five biggest cd4bs clusters
    print_stuff(annotations[cluster])

# bridging region epitope bnAb VH gene usage
def bridging_bool(q):
    if (annotations[cluster]['v_gene']).split('*')[0] not in bridging_genes:
        return False
    if annotations[cluster]['cdr3_length'] < 66:
        return False
    return True

bridging_clusters = [cluster for cluster in biggest_clusters if bridging_bool(cluster)]
print '\x1b[1;32;40m' + '  found %d bridging region clusters, printing the largest hits' % len(bridging_clusters)
sorted_bridging = sorted(bridging_clusters, key=lambda q: len(annotations[q]['unique_ids']), reverse=True)
for cluster in sorted_bridging[0:10]: # five biggest cd4bs clusters
    print_stuff(annotations[cluster])

# mper region epitope bnAb VH gene usage
def mper_bool(q):
    if (annotations[cluster]['v_gene']).split('*')[0] not in mper_genes:
        return False
    if annotations[cluster]['cdr3_length'] < 60:
        return False
    return True

mper_clusters = [cluster for cluster in biggest_clusters if mper_bool(cluster)]
print '\x1b[1;32;40m' + '  found %d mper region clusters, printing the largest hits' % len(mper_clusters)
sorted_mper = sorted(mper_clusters, key=lambda q: len(annotations[q]['unique_ids']), reverse=True)
for cluster in sorted_mper[0:10]: # five biggest cd4bs clusters
    print_stuff(annotations[cluster])

sys.exit()
