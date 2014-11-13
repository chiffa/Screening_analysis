__author__ = 'ank'

from csv import reader
from os import path
import numpy as np
from matplotlib import pyplot as plt
from chiffatools import hmm
from chiffatools. Linalg_routines import rm_nans, show_matrix_with_names, hierchical_clustering
from scipy.stats import ttest_ind
from itertools import combinations
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist

# intra-chromosome v.s. interchromosome variance?
# normalized within 50% lowest of the variance within a single chromosome?

pth = 'U:\\ank\\2014\\BRU_GBO\\4th_gen'
fle = 'mmc2-karyotypes.csv'

#####################################################
# Reads the source file for the amplification markers
#####################################################
with open(path.join(pth, fle)) as src:
    rdr = reader(src)
    header = rdr.next()
    selector = [1] + range(4, len(header))
    header = np.array(header)[selector]
    locuses = []
    for row in rdr:
        row = [ 'nan' if elt == 'NA' else elt for elt in row]
        locuses.append(np.genfromtxt(np.array(row))[selector].astype(np.float64))

################################
# Recovers the chromosome limits
################################
locuses = np.array(locuses)
chr_arr = locuses[:, 0]
broken_table = []
chroms = range(int(np.min(chr_arr)), int(np.max(chr_arr))+1)
print 'chroms:', chroms
for i in chroms:
    broken_table.append(chr_arr == i)

broken_table = np.array(broken_table)
chromosome_tag = np.repeat(locuses[:, 0].reshape((1, locuses.shape[0])), 200, axis=0)

#########################
# Defines the parsing HMM
#########################
transition_probs  = np.ones((3, 3)) * 0.01
np.fill_diagonal(transition_probs, 0.98)
initial_dist = np.array([[0.33, 0.34, 0.33]])
emission_probs  = np.ones((3, 3)) * 0.05
np.fill_diagonal(emission_probs, 0.9)
parsing_hmm = hmm.HMM(transition_probs, emission_probs)

def t_test_matrix(lane):

    current_lane = locuses[:, lane]
    dst = broken_table.shape[0]
    p_vals = np.empty((dst, dst))
    p_vals.fill(np.NaN)

    for i, j in combinations(range(0, dst), 2):
        _, p_val = ttest_ind(rm_nans(current_lane[broken_table[i, :]]), rm_nans(current_lane[broken_table[j, :]]), False)
        p_vals[i, j] = p_val
        # median = 1-np.abs(np.median(rm_nans(current_lane[broken_table[i, :]])) -
        #                 np.median(rm_nans(current_lane[broken_table[j, :]]))) / np.average(np.median(rm_nans(current_lane[broken_table[i, :]])) -
        #                                                                         np.median(rm_nans(current_lane[broken_table[j, :]])))
        # p_vals[i, j] = np.sqrt(p_val*median)
    return p_vals

# TODO: possible resolutions:
#       - mask out local amplification with HMM.
#       -

def compute_karyotype(lane, plotting=False, threshold = 0.33):

    def support_function(x, y):
        if x == 0 and y == 0:
            return 0
        if x == -1 and y <= 0:
            return -1
        if x >= 0 and y == 1:
            return 1
        if x == -1 and y == 1:
            return 0

    def plot_classification():

        classification_tag = np.repeat(parsed.reshape((1, parsed.shape[0])), 100, axis=0)
        ax1 = plt.subplot(311)
        plt.imshow(chromosome_tag, interpolation='nearest', cmap='spectral')
        plt.imshow(classification_tag, interpolation='nearest', cmap='coolwarm')
        plt.setp(ax1.get_xticklabels(), fontsize=6)
        ax2 = plt.subplot(312, sharex=ax1)
        plt.plot(current_lane, 'k.', label = '%s %s %s'%(lane,
                                                         '{0:.2f}'.format(np.mean(rm_nans(current_lane))),
                                                         '{0:.2f}'.format(np.std(rm_nans(current_lane)))))
        plt.legend(prop={'size':10})
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax3 = plt.subplot(313, sharex=ax1)
        plt.imshow(chromosome_tag, interpolation='nearest', cmap='spectral')
        plt.imshow(re_classification_tag, interpolation='nearest', cmap='coolwarm')
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.show()

    current_lane = locuses[:, lane]
    current_lane = current_lane - np.mean(rm_nans(current_lane))
    current_lane = current_lane / np.std(rm_nans(current_lane))*0.3
    binarized = (current_lane > threshold).astype(np.int16) - (current_lane < -threshold) + 1
    parsed = np.array(hmm.viterbi(parsing_hmm, initial_dist, binarized))

    t_mat = t_test_matrix(lane)
    t_mat[np.isnan(t_mat)] = 0
    t_mat = t_mat + t_mat.T
    np.fill_diagonal(t_mat, 1)
    ct_mat = t_mat.copy()
    ct_mat[t_mat < 0.01] = 0.01
    ct_mat = 1 - ct_mat

    Y = sch.linkage(ct_mat, method='centroid')
    clust_alloc = sch.fcluster(Y, 0.95, criterion='distance')

    averages = []
    for i in range(0, 24):
        averages.append(np.median(rm_nans(current_lane[broken_table[i]])))

    accumulator = [[] for _ in range(0, max(clust_alloc)+1)]
    for loc, item in enumerate(averages):
        accumulator[clust_alloc[loc]].append(item)
    accumulator = np.array([ np.average(np.array(_list)) for _list in accumulator][1:])

    splitting_matrix = np.repeat(accumulator.reshape((1, accumulator.shape[0])),
                                    accumulator.shape[0], axis = 0)
    splitting_matrix = np.abs(splitting_matrix - splitting_matrix.T)

    # plt.imshow(splitting_matrix, interpolation='nearest', cmap='coolwarm')
    # plt.show()

    Y_2 = sch.linkage(splitting_matrix, method='centroid')
    clust_alloc_2 = sch.fcluster(Y_2, 3, criterion='maxclust')

    accumulator_2 = [[] for _ in range(0, max(clust_alloc_2)+1)]
    for loc, item in enumerate(accumulator):
        accumulator_2[clust_alloc_2[loc]].append(item)
    accumulator_2 = np.array([ np.average(np.array(_list)) for _list in accumulator_2][1:])

    sorter_l = np.argsort(accumulator_2)
    sorter = dict((pos, i) for i, pos in enumerate(sorter_l))

    re_chromosome_pad = np.repeat(locuses[:, 0].reshape((1, locuses.shape[0])), 100, axis=0)
    re_classification_tag = np.zeros(re_chromosome_pad.shape)

    for i in range(0, len(clust_alloc)):
        # re_classification_tag[re_chromosome_pad == i+1] = averages[i]
        # re_classification_tag[re_chromosome_pad == i+1] = accumulator[ clust_alloc[i]-1]
        re_classification_tag[re_chromosome_pad == i+1] = sorter[clust_alloc_2[clust_alloc[i]-1]-1]

    collector = []
    for i in chroms:
        lw = np.percentile(parsed[broken_table[i-1, :]]-1, 25)
        hr = np.percentile(parsed[broken_table[i-1, :]]-1, 75)
        collector.append(support_function(lw, hr))

    if plotting:
        plot_classification()

    collector = np.array(collector)
    return collector


def compute_all_karyotypes():

    def plot():
        plt.imshow(chromlist, interpolation='nearest', cmap='coolwarm')
        plt.show()

    chromlist = []
    for i in range(1, locuses.shape[1]):
        chromlist.append(compute_karyotype(i, plotting=True, threshold=0.35))
    chromlist = np.array(chromlist).astype(np.float64)
    return chromlist, header[1:locuses.shape[1]]


if __name__ == "__main__":
    # compute_karyotype(6, True, 0.35)
    print compute_all_karyotypes()