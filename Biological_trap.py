__author__ = 'ank'

import numpy as np
from csv import reader
from os import path
from matplotlib import pyplot as plt
from math import sqrt
from Support import pretty_gradual_plot, get_resistant_susceptible
import matplotlib as mlb
from itertools import permutations
from Linalg_routines import hierchical_clustering, show_matrix_with_names

pth = 'U:\\ank\\2014\\BRU_GBO\\4th_gen'
fle = 'gb-breast_cancer.tsv'


headers = ['cellline', 'compound', 'drug_plate_id', 'T0_plate_id', 'background_od1', 'background_od2',
           'od0.1', 'od0.2', 'od0.3', 'od1.1', 'od1.2', 'od1.3', 'od2.1', 'od2.2', 'od2.3', 'od3.1',
           'od3.2', 'od3.3', 'od4.1', 'od4.2', 'od4.3', 'od5.1', 'od5.2', 'od5.3', 'od6.1', 'od6.2',
           'od6.3', 'od7.1', 'od7.2', 'od7.3', 'od8.1', 'od8.2', 'od8.3', 'od9.1', 'od9.2', 'od9.3',
           'T0_background_od1', 'T0_background_od2', 'T0_median_od', 'c1', 'c2', 'c3', 'c4', 'c5',
           'c6', 'c7', 'c8', 'c9', 'units']


def index(myset):
    return dict((elt, i) for i, elt in enumerate(myset))

def broadcast(subline):
    if len(subline) !=30:
        print subline
        raise Exception('wrong number of items in subline')
    else:
        arr = np.array(subline)
        return arr.reshape((10, 3))

def remove_nan(array):
    array = array.flatten()
    return array[np.logical_not(np.isnan(array))]

def make_comparator(percentile_5_range):
    st = sqrt(2)

    def compare(val1, val2):
        if val1-val2 > st*percentile_5_range:
            return 1
        if val1-val2 < -st*percentile_5_range:
            return -1
        else:
            return 0

    return compare


with open(path.join(pth, fle)) as src:
    cells = []
    drugs = []
    rdr = reader(src, dialect='excel-tab')
    header = rdr.next()
    for row in rdr:
        cells.append(row[0])
        drugs.append(row[1])

cell_idx = index(set(cells))
drug_idx = index(set(drugs))

cell_idx_rv = dict([(value, key) for key, value in cell_idx.iteritems()])
drug_idx_rv = dict([(value, key) for key, value in drug_idx.iteritems()])

cellno = len(cell_idx)
drugno = len(drug_idx)

# print cell_idx
# print drug_idx

# print cellno, drugno

storage = np.empty((cellno, drugno, 10, 3))
storage.fill(np.NaN)

background = np.empty((cellno, drugno, 4))
background.fill(np.NaN)

concentrations = np.empty((cellno, drugno, 10))
concentrations.fill(np.NaN)

with open(path.join(pth, fle)) as src:
    rdr = reader(src, dialect='excel-tab')
    test_array = rdr.next()
    broadcast(test_array[6:36])
    for row in rdr:
        cell_no = cell_idx[row[0]]
        drug_no = drug_idx[row[1]]
        storage[cell_no, drug_no, :, :] = broadcast(row[6:36])
        background[cell_no, drug_no, :] = np.array([row[i] for i in [4, 5, 36, 37]])
        concentrations[cell_no, drug_no, :] = np.array([0]+row[39:48])

# print storage
# print background
# print concentrations
# print np.std(remove_nan(background))
median = np.percentile(remove_nan(background), 50)
sensible_min = np.percentile(remove_nan(background), 2)
sensible_max = np.percentile(remove_nan(background), 98)

compare = make_comparator(max(sensible_max - median, median - sensible_min))

mlb.rcParams['figure.figsize'] = (20,15)

drug = 'Vorinostat'
example = storage[:, drug_idx[drug], :, :]
example_concs = concentrations[:, drug_idx[drug], :]

# pretty_gradual_plot(example, example_concs, cell_idx_rv, drug, blank_line=sensible_max)
# print get_resistant_susceptible(example, drug, blank_line=sensible_max)


foldlist = []
for index, drug in sorted(drug_idx_rv.items()):
    drug_action = storage[:, drug_idx[drug], :, :]
    foldlist.append(get_resistant_susceptible(drug_action, drug, blank_line=sensible_max))
flds = np.array(foldlist)



plt.imshow(flds, interpolation='nearest', cmap='coolwarm')
plt.yticks(range(0, len(drug_idx)), [drug for i, drug in sorted(drug_idx_rv.items())], rotation='horizontal')
plt.xticks(range(0, len(cell_idx)), [cell for i, cell in sorted(cell_idx_rv.items())], rotation='vertical')
plt.colorbar()
plt.show()

from scipy.stats import itemfreq
fulter = np.all(np.logical_not(np.isnan(storage)), axis=(2,3))
cnt = np.sum(fulter, axis = 1)
for i, cell in sorted(cell_idx_rv.items()):
    ifrq = itemfreq(remove_nan(flds[:,i])+1)
    collector = [0, 0, 0]
    for val, count in ifrq:
        collector[int(val)]=int(count/float(cnt[i])*100)
    if len(cell)<6:
        pad = ''.join([' ']*(6-len(cell)))
    else:
        pad = ''
    print cell+pad, ':', '\t', 'kill: %s %% \t weak: %s %% \t strong: %s %% '% (collector[0], collector[1], collector[2])


# todo: plot the drugs - cell combination
#
# count = np.sum(fulter, axis = 0).astype(np.float)
# print count.shape
# p_drug = np.sum(flds<0.01, axis=1)/count
# collector_array = np.zeros((p_drug.shape[0], p_drug.shape[0]))
# for i, j in permutations(range(0, p_drug.shape[0]), 2):
#     # collector_array[i, j] = np.sum(np.logical_or(flds[i, :]<0.01, flds[j, :]<0.01))/float(max(count[i], count[j]))/np.sqrt(p_drug[i]*p_drug[j])
#     collector_array[i, j] = max(1.5, float(np.sum(np.logical_xor(flds[i, :]<0.01, flds[j, :]<0.01))-np.sum(np.logical_xor(fulter[:, i], fulter[:, j]))))
#
# collector_array = 1./(collector_array+np.identity(collector_array.shape[0]))
# drugnames = [drug for i, drug in sorted(drug_idx_rv.items())]
# index = hierchical_clustering(collector_array, drugnames)
#
# drugnames = np.array(drugnames)[index].tolist()
# collector_array = collector_array[:, index]
# collector_array = collector_array[index, :]
#
# show_matrix_with_names(1./collector_array, drugnames, drugnames)