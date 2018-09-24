# Oct 2017
# A rewrite of the original gene usage plotting scripts

##### Pseudocode ######
# 1. read in .freq files (start with hv only)
# 2. depending on user-specified gene option, store relevant detail, including freq
# 3. count proportion of either unique or total TCRs seen in sample
# 4. plots distribution using violinplot..? or boxplot..?

##### Py packages #####
from __future__ import division
from os import listdir
from os.path import isfile, join
from collections import defaultdict
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
import urllib2 as url
import numpy as np
import gzip
import pylab
import datetime
import sys

##### Fx blocks #####

def getVals(chain, file_dir, dcr_id):

    files = [f for f in listdir(file_dir) if isfile(join(file_dir, f)) if chain in f]
    
    comb_gene_prop = defaultdict(list)

    for f in files:
        with gzip.open(file_dir+f, 'r') as inf:
            lines = inf.read().splitlines()
            freq = [int(l.split(', ')[5]) for l in lines]

            gene_freq = defaultdict(list)

            for l in lines:
                spl = l.split(', ')
                gene_freq[spl[int(dcr_id)]].append(int(spl[5]))

            for g, freq_lst in gene_freq.iteritems():
                # total
                #prop = sum(freq_lst)/sum(freq)
                # or unique
                prop = len(freq_lst)/len(freq)
                comb_gene_prop[g].append(prop)

    ref_tag_file, chr_order = getGeneNameAndOrder(gene)
    num = 0
    ref_tag = {}
    
    for line in ref_tag_file:
        name = line.split('|')[1]
        ref_tag[num] = name
        num += 1
        
    gene_vals = {}
    labels = []
    y = [[0]]*len(chr_order)

    for idx, order in enumerate(chr_order):
        for g1, name in ref_tag.iteritems():
            if order == g1:
                labels.append(name)
                for g2, prop in comb_gene_prop.iteritems():
                    if str(g1) == g2:
                        y[idx] = prop

    return y, labels

def getGeneNameAndOrder(gene):

    if gene == 'TRAV':
        add_to_trav = 'https://raw.githubusercontent.com/innate2adaptive/Decombinator-Tags-FASTAs/master/human_extended_TRAV.tags'
        ref_tag_file = url.urlopen(add_to_trav)

        # missing : 8-4
        trav_chr_order = [0, 1, 13, 24, 32, 35, 36, 38, 2, 3, 39, 40, 6, 4, 7, 8, 43, 5, 41, 9, 10, 12, 14, 15, 16, 17, 18, 19, 20, 22, 23, 25, 21, 26, 27, 28, 29, 30, 31, 33, 34]

        return ref_tag_file, trav_chr_order

    if gene == 'TRAJ':
        add_to_traj = 'https://raw.githubusercontent.com/innate2adaptive/Decombinator-Tags-FASTAs/master/human_extended_TRAJ.tags'
        ref_tag_file = url.urlopen(add_to_traj)

        traj_chr_order = [60, 59, 58, 57, 45, 44, 56, 43, 42, 41, 55, 40, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 27, 26, 25, 24, 54, 23, 22, 21, 20, 19, 17, 16, 15, 14, 53, 13, 12, 11, 10, 9, 52, 8, 7, 6, 5, 4, 3, 2, 1, 0, 49, 48, 47, 46, 39, 28, 18, 51, 50]

        return ref_tag_file, traj_chr_order

    if gene == 'TRBV':
        add_to_trbv = 'https://raw.githubusercontent.com/innate2adaptive/Decombinator-Tags-FASTAs/master/human_extended_TRBV.tags'
        ref_tag_file = url.urlopen(add_to_trbv)

        # missing 6-3, 3-2, 6-2, 12-3
        trbv_chr_order = [51, 14, 21, 23, 26, 31, 50, 24, 25, 37, 59, 57, 32, 38, 60, 47, 44, 0, 3, 52, 1, 4, 53, 33, 39, 27, 34, 58, 28, 49, 40, 29, 35, 41, 48, 36, 42, 30, 43, 8, 2, 5, 6, 7, 9, 10, 11, 45, 12, 13, 15, 54, 46, 16, 17, 61, 56, 62, 18, 19, 20, 22]

        return ref_tag_file, trbv_chr_order

    if gene == 'TRBJ':
        add_to_trbj = 'https://raw.githubusercontent.com/innate2adaptive/Decombinator-Tags-FASTAs/master/human_extended_TRBJ.tags'
        ref_tag_file = url.urlopen(add_to_trbj)

        trbj_chr_order = [0, 1, 2, 3, 4, 5, 6, 7, 13, 8, 9, 10, 11, 12]

        return ref_tag_file, trbj_chr_order
    
##### Start #####

chain = str(sys.argv[1])
gene = str(sys.argv[2])
dcr_id = str(sys.argv[3])
#hv_dir = '/Volumes/BF_MI_1/cvid-analysis/hv_freq_chainchud/'
hv_dir = '/Volumes/BF_MI_1/cvid-analysis/hv_freq/'
id_dir = '/Volumes/BF_MI_1/cvid-analysis/cvid_freq_notempus/'

y1, labels1 = getVals(chain, hv_dir, dcr_id)
y2, labels2 = getVals(chain, id_dir, dcr_id)

# create a figure instance
fig = plt.figure(1, figsize = (6, 8))

num = 1
for val in y2:
    y = np.repeat(num, len(val))
    plt.scatter(val, y, color='green', s=10)
    num += 1

# box plot
# create an axes instance
ax = fig.add_subplot(111)

# create the boxplot
bp1 = ax.boxplot(y1, whis='range', vert=0)
pylab.setp(bp1['whiskers'], color = 'black')
pylab.setp(bp1['fliers'], marker = 'None')
pylab.setp(bp1['medians'], color = 'black')
pylab.setp(bp1['boxes'], color = 'black')

plt.yticks(range(1, len(y1)+1, 1), labels1, rotation='horizontal', fontsize=10)
plt.xticks(rotation=45)
plt.gca().set_xlim(left=-0.01)
plt.xlabel('Proportion (total unique)')
plt.tight_layout()
#plt.show()
today_date = datetime.datetime.now()
plt.savefig("%s-%s-%s_" %(today_date.year, today_date.month, today_date.day)+str(sys.argv[2])+"-UniqueProp-HVvsID.png", dpi=600)
