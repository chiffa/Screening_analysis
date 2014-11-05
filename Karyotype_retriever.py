__author__ = 'ank'

from csv import reader
from os import path
import numpy as np
from matplotlib import pyplot as plt
from chiffatools import hmm
# from chiffatools.Linalg_routines import s
from scipy.stats import gaussian_kde

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


def compute_karyotype(lane, plotting=False, threshold = 0.33):

    def support_function(x, y):
        if x == 0 and y == 0:
            return 0
        if x == -1 and y <= 0:
            return -1
        if x >= 0 and y == 1:
            return 1
        if x == -1 and y == 1:
            return np.NaN

    def plot_classification():

        classification_tag = np.repeat(parsed.reshape((1, parsed.shape[0])), 100, axis=0)
        ax1 = plt.subplot(311)
        plt.imshow(chromosome_tag, interpolation='nearest', cmap='spectral')
        plt.imshow(classification_tag, interpolation='nearest', cmap='coolwarm')
        plt.setp(ax1.get_xticklabels(), fontsize=6)
        ax2 = plt.subplot(312, sharex=ax1)
        plt.plot(current_lane, 'k.')
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.show()

    current_lane = locuses[:, lane]
    binarized = (current_lane > threshold).astype(np.int16) - (current_lane < -threshold) + 1
    parsed = np.array(hmm.viterbi(parsing_hmm, initial_dist, binarized))

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
        chromlist.append(compute_karyotype(i, plotting=False, threshold=0.35))
    chromlist = np.array(chromlist).astype(np.float64)
    print chromlist
    plot()
    return chromlist

if __name__ == "__main__":
    compute_all_karyotypes()