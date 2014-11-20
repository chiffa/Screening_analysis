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
from scipy.ndimage.filters import gaussian_filter1d as smooth_signal
import warnings

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
transition_probs  = np.ones((3, 3)) * 0.001
np.fill_diagonal(transition_probs, 0.998)
initial_dist = np.array([[0.33, 0.34, 0.33]])
emission_probs  = np.ones((3, 3)) * 0.1
np.fill_diagonal(emission_probs, 0.8)
parsing_hmm = hmm.HMM(transition_probs, emission_probs)
parsing_hmm2 = hmm.HMM(transition_probs, emission_probs)

warnings.catch_warnings()
warnings.simplefilter("ignore")

def t_test_matrix(lane):

    current_lane = locuses[:, lane]
    dst = broken_table.shape[0]
    p_vals = np.empty((dst, dst))
    p_vals.fill(np.NaN)

    for i, j in combinations(range(0, dst), 2):
        _, p_val = ttest_ind(rm_nans(current_lane[broken_table[i, :]]), rm_nans(current_lane[broken_table[j, :]]), False)
        p_vals[i, j] = p_val
    return p_vals

def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


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
        classification_tag2 = np.repeat(parsed2.reshape((1, parsed.shape[0])), 100, axis=0)
        ax1 = plt.subplot(311)
        plt.imshow(chromosome_tag, interpolation='nearest', cmap='spectral')
        plt.imshow(classification_tag, interpolation='nearest', cmap='coolwarm')
        plt.setp(ax1.get_xticklabels(), fontsize=6)
        ax2 = plt.subplot(312, sharex=ax1)
        plt.plot(current_lane, 'k.', label = '%s %s %s'%(lane,
                                                         '{0:.2f}'.format(np.mean(rm_nans(current_lane))),
                                                         '{0:.2f}'.format(np.std(rm_nans(current_lane)))))
        plt.plot(gauss_convolve, 'r', lw=2)
        plt.plot(gauss_convolve+rolling_std, 'g', lw=2)
        plt.plot(gauss_convolve-rolling_std, 'g', lw=2)
        plt.legend(prop={'size':10})
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax3 = plt.subplot(313, sharex=ax1)
        # plt.plot(total_corr, 'k.')
        # plt.plot(gauss_convolve, 'k.')
        # plt.plot(rolling_std, 'g')
        plt.imshow(chromosome_tag, interpolation='nearest', cmap='spectral')
        plt.imshow(classification_tag2, interpolation='nearest', cmap='coolwarm')
        # plt.setp(ax3.get_xticklabels(), visible=False)
        plt.show()

    current_lane = locuses[:, lane]
    current_lane = current_lane - np.mean(rm_nans(current_lane))
    current_lane = current_lane / np.std(rm_nans(current_lane))*0.3
    binarized = (current_lane > threshold).astype(np.int16) - (current_lane < -threshold) + 1
    parsed = np.array(hmm.viterbi(parsing_hmm, initial_dist, binarized))

    partial_corr = np.correlate(rm_nans(current_lane), rm_nans(current_lane), mode='full')
    total_corr = partial_corr[partial_corr.size/2:]

    gauss_convolve = np.empty(current_lane.shape)
    gauss_convolve.fill(np.NaN)
    gauss_convolve[np.logical_not(np.isnan(current_lane))] = smooth_signal(rm_nans(current_lane), 10, order=0, mode='mirror')

    rolling_std = np.empty(current_lane.shape)
    rolling_std.fill(np.NaN)
    payload = np.std(rolling_window(rm_nans(current_lane), 10), 1)
    c1, c2 = (np.sum(np.isnan(current_lane[:5])), np.sum(np.isnan(current_lane[-4:])))
    rolling_std[5:-4][np.logical_not(np.isnan(current_lane))[5:-4]] = np.lib.pad(payload,
                                                                                 (c1, c2),
                                                                                 'constant',
                                                                                 constant_values=(np.NaN, np.NaN))

    # nns_rolling_std = rm_nans(rolling_std)
    # print(np.percentile(nns_rolling_std, 50),
    #       np.percentile(nns_rolling_std, 75),
    #       np.percentile(nns_rolling_std, 80),
    #       np.percentile(nns_rolling_std, 90),
    #       np.percentile(nns_rolling_std, 95))

    prct_threshold = np.percentile(rm_nans(rolling_std), 75)
    print prct_threshold
    binarized2 = (current_lane > prct_threshold).astype(np.int16) - (current_lane < -prct_threshold) + 1
    parsed2 = np.array(hmm.viterbi(parsing_hmm2, initial_dist, binarized2))

    # print broken_table.shape[0]
    # for i in range(0, broken_table.shape[0]):
    #     current_chromosome = rm_nans(current_lane[broken_table[i, :]])
    #     print i+1, np.mean(current_chromosome), np.std(current_chromosome)
    #     print '\t', np.median(current_chromosome), (np.percentile(current_chromosome, 85)-np.percentile(current_chromosome, 15))/2.0
    #


    t_mat = t_test_matrix(lane)
    t_mat[np.isnan(t_mat)] = 0
    t_mat = t_mat + t_mat.T
    np.fill_diagonal(t_mat, 1)
    ct_mat = t_mat.copy()
    ct_mat[t_mat < 0.01] = 0.01
    ct_mat = 1 - ct_mat

    # TODO: idea: keep the HMM model for the regression, but determine bounds by normalizing
    #   On the per-chromosome basis: compute the mean and std
    #   Eliminate the chromosomes with highest std (possibly contain inner amplifications)
    #   determine the bounds of the HMM model.
    #   Run the HMM rounds iteratively,
    #       On each iteration, remove losses/gains from the previous model.
    #       End up with discretized levels.
    #   Try to gess most likely integer levels of chromosomes/locus gains.

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