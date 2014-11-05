__author__ = 'ank'

from csv import reader
from os import path
import numpy as np
from matplotlib import pyplot as plt
from chiffatools import hmm


pth = 'U:\\ank\\2014\\BRU_GBO\\4th_gen'
fle = 'mmc2-karyotypes.csv'


with open(path.join(pth, fle)) as src:
    rdr = reader(src)
    header = rdr.next()
    selector = [1] + range(4, len(header))
    header = np.array(header)[selector]
    locuses = []
    for row in rdr:
        row = [ 'nan' if elt == 'NA' else elt for elt in row]
        locuses.append(np.genfromtxt(np.array(row))[selector].astype(np.float64))


locuses = np.array(locuses)
chromosome_tag = np.repeat(locuses[:, 0].reshape((1, locuses.shape[0])), 200, axis=0)

binarized = (locuses[:, 2] > 0.33).astype(np.int16) - (locuses[:, 2] < -0.33) + 1
transition_probs  = np.ones((3, 3)) * 0.01
np.fill_diagonal(transition_probs, 0.98)
initial_dist = np.array([[0.33, 0.34, 0.33]])
emission_probs  = np.ones((3, 3)) * 0.05
np.fill_diagonal(emission_probs, 0.9)
parsing_hmm = hmm.HMM(transition_probs, emission_probs)
parsed = np.array(hmm.viterbi(parsing_hmm, initial_dist, binarized))
parsed = parsed.reshape((1, parsed.shape[0]))

classification_tag = np.repeat(parsed, 100, axis=0)
print parsed.shape
# chromosome_tag = np.column_stack(chromosome_tag, classification_tag)

ax1 = plt.subplot(211)
plt.imshow(chromosome_tag, interpolation='nearest', cmap='spectral')
plt.imshow(classification_tag, interpolation='nearest', cmap='coolwarm')
plt.setp(ax1.get_xticklabels(), fontsize=6)
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(locuses[:, 2], 'k.')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.show()