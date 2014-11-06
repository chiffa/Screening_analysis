__author__ = 'ank'

import numpy as np
from csv import reader
from os import path
from matplotlib import pyplot as plt
from math import sqrt
from Support import pretty_gradual_plot, get_resistant_susceptible, pretty_get_resistant_susceptible
import matplotlib as mlb
from itertools import product
from chiffatools.Linalg_routines import hierchical_clustering, show_matrix_with_names
from scipy.stats import gaussian_kde
from matplotlib.pyplot import Rectangle
from Karyotype_retriever import compute_all_karyotypes
from sklearn import linear_model
from scipy.stats import spearmanr
from chiffatools.Linalg_routines import show_matrix_with_names

locus, cell_line_name = compute_all_karyotypes()

pth = 'U:\\ank\\2014\\BRU_GBO\\4th_gen'
fle = 'gb-breast_cancer.tsv'


headers = ['cellline', 'compound', 'drug_plate_id', 'T0_plate_id', 'background_od1', 'background_od2',
           'od0.1', 'od0.2', 'od0.3', 'od1.1', 'od1.2', 'od1.3', 'od2.1', 'od2.2', 'od2.3', 'od3.1',
           'od3.2', 'od3.3', 'od4.1', 'od4.2', 'od4.3', 'od5.1', 'od5.2', 'od5.3', 'od6.1', 'od6.2',
           'od6.3', 'od7.1', 'od7.2', 'od7.3', 'od8.1', 'od8.2', 'od8.3', 'od9.1', 'od9.2', 'od9.3',
           'T0_background_od1', 'T0_background_od2', 'T0_median_od', 'c1', 'c2', 'c3', 'c4', 'c5',
           'c6', 'c7', 'c8', 'c9', 'units']


def index_f(myset):
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

cell_idx = index_f(set(cells))
drug_idx = index_f(set(drugs))

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


# print storage.shape
# print background
# print concentrations
# print np.std(remove_nan(background))
median = np.percentile(remove_nan(background), 50)
sensible_min = np.percentile(remove_nan(background), 2)
sensible_max = np.percentile(remove_nan(background), 98)

compare = make_comparator(max(sensible_max - median, median - sensible_min))

mlb.rcParams['figure.figsize'] = (20, 15)

# drug = '17-AAG'
# example = storage[:, drug_idx[drug], :, :]
# example_concs = concentrations[:, drug_idx[drug], :]
#
# pretty_gradual_plot(example, example_concs, cell_idx_rv, drug, blank_line=sensible_max)
# print get_resistant_susceptible(example, drug, blank_line=sensible_max)

foldlist = []
for index, drug in sorted(drug_idx_rv.items()):
    drug_action = storage[:, drug_idx[drug], :, :]
    foldlist.append(get_resistant_susceptible(drug_action, drug, blank_line=sensible_max))
    # if 'Tykerb' in drug:
    #     pretty_gradual_plot(storage[:, index, :, :], concentrations[:, index, :], cell_idx_rv, drug, blank_line=sensible_max)
    #     pretty_get_resistant_susceptible(drug_action, cell_idx_rv, drug,
    #             blank_line=sensible_max, concentrations=concentrations[:, drug_idx[drug], :])
flds = np.array(foldlist)


from scipy.stats import itemfreq
fulter = np.all(np.logical_not(np.isnan(storage)), axis=(2,3))
cnt = np.sum(fulter, axis = 1)
for i, cell in sorted(cell_idx_rv.items()):
    ifrq = itemfreq((remove_nan(flds[:,i])+1)*2)
    collector = [0, 0, 0, 0, 0]
    for val, count in ifrq:
        collector[int(val)]=int(count/float(cnt[i])*100)
    if len(cell)<10:
        pad = ''.join([' ']*(10-len(cell)))
    else:
        pad = ''
    print cell+pad, ':', '\t', 'kill: %s %% \t susceptible: %s %% \t partially resistant: %s %% \t resistant: %s %%'% (collector[0], collector[1], collector[3], collector[4])


cm = plt.cm.get_cmap('coolwarm')
red = Rectangle((0, 0), 1, 1, fc=cm(1.0))
light_red = Rectangle((0, 0), 1, 1, fc=cm(0.75))
grey = Rectangle((0, 0), 1, 1, fc=cm(0.5))
light_blue = Rectangle((0, 0), 1, 1, fc=cm(0.25))
blue = Rectangle((0, 0), 1, 1, fc=cm(0.0))


def characterize_array(matrix, title, analysis_type=0):
    experiment_mask = np.logical_not(np.isnan(matrix))
    collector_array = np.zeros((matrix.shape[0], matrix.shape[0]))
    for i, j in product(range(0, matrix.shape[0]), repeat=2):
        common_ground = np.logical_and(experiment_mask[i, :], experiment_mask[j, :])
        action1 = float(np.sum(np.logical_and(matrix[i, :]<-0.1, common_ground)))/float(np.sum(common_ground))
        action2 = float(np.sum(np.logical_and(matrix[j, :]<-0.1, common_ground)))/float(np.sum(common_ground))
        if analysis_type:
            complementary_action = float(np.sum(np.logical_and(np.logical_or(matrix[i, :]<-0.1, matrix[j, :]<-0.1), common_ground)))/float(np.sum(common_ground))
            collector_array[i, j] = complementary_action
        else:
            complementary_action = float(np.sum(np.logical_and(np.logical_xor(matrix[i, :]<-0.1, matrix[j, :]<-0.1), common_ground)))/float(np.sum(common_ground))
            collector_array[i, j] = (complementary_action)/(1 - (action1 * action2 + (1 - action1) * (1 - action2)))

    density = gaussian_kde(collector_array.flatten())
    xs = np.linspace(collector_array.min(), collector_array.max(), 200)
    plt.hist(collector_array.flatten(), 100, histtype='step', normed=True, label=title)
    plt.plot(xs, density(xs), label=title)

    return collector_array


def random_permute(matrix, axis):

    def na_stable_permutation(D1_array):
        msk = np.logical_not(np.isnan(D1_array))
        load = np.random.permutation(D1_array[msk])
        new_arr = np.empty(D1_array.shape)
        new_arr.fill(np.NaN)
        new_arr[msk] = load
        print new_arr
        print D1_array
        return new_arr

    return np.apply_along_axis(na_stable_permutation, axis=axis, arr=matrix)


def plot_action_map(fields, Drug1=None, Drug2=None):

    if Drug1 is None and Drug2 is None:
        # min / max
        # percent killed at max concentration (drug)
        # percent resistant or weakly resistant
        # percent killed or susceptible
        plt.imshow(fields, interpolation='nearest', cmap='coolwarm')
        plt.tick_params(axis='both',  labelsize=8)
        plt.yticks(range(0, len(drug_idx)), [drug for i, drug in sorted(drug_idx_rv.items())], rotation='horizontal', )
        plt.xticks(range(0, len(cell_idx)), [cell for i, cell in sorted(cell_idx_rv.items())], rotation='vertical')
        plt.legend((red, light_red, grey, light_blue, blue),
                   ('Resistant', 'Partially Resistant', 'Weakly resistant', 'Susceptible', 'Killed'),
                   bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.,)
        plt.show()

    if Drug1 is not None and Drug2 is not None:
        if type(Drug1) is str:
            idx1 = drug_idx[Drug1]
            name1 = Drug1
        if type(Drug1) is int:
            idx1 = Drug1
            name1 = drug_idx_rv[Drug1]
        if type(Drug2) is str:
            idx2 = drug_idx[Drug2]
            name2 = Drug2
        if type(Drug2) is int:
            idx2 = Drug2
            name2 = drug_idx_rv[Drug2]

        print idx1, name1, idx2, name2
        matrix = flds
        experiment_mask = np.logical_not(np.isnan(matrix))
        common_ground = np.logical_and(experiment_mask[idx1, :], experiment_mask[idx2, :])
        action1 = float(np.sum(np.logical_and(matrix[idx1, :]<-0.1, common_ground)))/float(np.sum(common_ground))
        action2 = float(np.sum(np.logical_and(matrix[idx2, :]<-0.1, common_ground)))/float(np.sum(common_ground))
        print action1, action2, float(np.sum(common_ground))
        print float(np.sum(np.logical_and(np.logical_or(matrix[idx1, :]<-0.1, matrix[idx2, :]<-0.1), common_ground)))
        complementary_action = float(np.sum(np.logical_and(np.logical_or(matrix[idx1, :]<-0.1, matrix[idx2, :]<-0.1), common_ground)))/float(np.sum(common_ground))
        print (complementary_action)
        complementary_action = float(np.sum(np.logical_and(np.logical_xor(matrix[idx1, :]<-0.1, matrix[idx2, :]<-0.1), common_ground)))/float(np.sum(common_ground))
        print (complementary_action)/(1 - (action1 * action2 + (1 - action1) * (1 - action2)))


        plt.imshow(fields[(idx1, idx2),:], interpolation='nearest', cmap='coolwarm')
        plt.tick_params(axis='both',  labelsize=8)
        plt.yticks(range(0, 2), [name1, name2], rotation='horizontal', )
        plt.xticks(range(0, len(cell_idx)), [cell for i, cell in sorted(cell_idx_rv.items())], rotation='vertical')
        plt.legend((red, light_red, grey, light_blue, blue),
                   ('Resistant', 'Partially Resistant', 'Weakly resistant', 'Susceptible', 'Killed'),
                   bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.,)
        plt.show()


def round_show(matrix):
    mat_source = matrix.copy()
    drugnames = [drug for i, drug in sorted(drug_idx_rv.items())]
    index = hierchical_clustering(matrix, drugnames)
    drugnames = np.array(drugnames)[index].tolist()
    matrix = matrix[:, index]
    matrix = matrix[index, :]
    matrix[matrix==101] = np.NaN
    # np.fill_diagonal(matrix, 1)
    # matrix = 1.0 / matrix
    show_matrix_with_names(matrix, drugnames, drugnames)

    matrix = mat_source
    matrix[matrix==101] = np.NaN
    np.fill_diagonal(matrix, 0)
    fltr = np.isnan(matrix)
    matrix[fltr] = 0
    crds = []
    while np.amax(matrix) > 0:
        coords = np.unravel_index(np.argmax(matrix), matrix.shape)
        crds.append(coords)
        matrix[coords] = 0
        matrix[coords[1], coords[0]] = 0

    for i, coords in enumerate(crds):
        idx1 = coords[0]
        name1 = drug_idx_rv[coords[0]]
        idx2 = coords[1]
        name2 = drug_idx_rv[coords[1]]

        ax = plt.subplot(len(crds), 1, i+1)
        ax.imshow(flds[(idx1, idx2),:], interpolation='nearest', cmap='coolwarm')
        ax.tick_params(axis='both',  labelsize=8)
        plt.yticks(range(0, 2), [name1, name2], rotation='horizontal', )
        if i==0:
            ax.legend((red, light_red, grey, light_blue, blue),
                   ('Resistant', 'Partially Resistant', 'Weakly resistant', 'Susceptible', 'Killed'),
                   bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.,)
        if i == len(crds)-1:
            plt.xticks(range(0, len(cell_idx)), [cell for i, cell in sorted(cell_idx_rv.items())], rotation='vertical')
    plt.show()
    # plot_action_map(flds, int(coords[0]), int(coords[1]))
    return index



# plot_action_map(flds)
# plot_action_map(flds, 'Rapamycin', 'Fascaplysin')

# per_flds = random_permute(flds, 1)
# plot_action_map(per_flds)

collector_array = characterize_array(flds, 'original', 1)
plt.clf()
# drugnames = [drug for i, drug in sorted(drug_idx_rv.items())]
# show_matrix_with_names(collector_array, drugnames, drugnames)
# round_show(collector_array)
# permcollector_array = characterize_array(per_flds, 'permuted')
# plt.legend()
# plt.show()
# plt.clf()

###############################################################################
###############################################################################
###############################################################################

print locus.shape
cell_line_name = cell_line_name.tolist()
c_dict = index_f(cell_line_name)
print cell_line_name
print cell_idx.keys()

ordered_locus = np.empty((24, cellno))
print ordered_locus.shape
ordered_locus.fill(np.NaN)

count = 0
for i, name in cell_idx_rv.iteritems():
    name_idx = c_dict.get(name, -1)
    if name_idx < 0:
        print '%s not found' % name
        count += 1
    else:
        ordered_locus[:, i] = locus[name_idx, :]

print count
print ordered_locus
# plt.imshow(ordered_locus, interpolation='nearest', cmap='coolwarm')
# plt.show()

# drug = '17-AAG'
# example = storage[:, drug_idx[drug], :, :]
# example_concs = concentrations[:, drug_idx[drug], :]


ax1 = plt.subplot(211)
plt.imshow(flds, interpolation='nearest', cmap='coolwarm')
plt.setp(ax1.get_xticklabels(), fontsize=6)
ax2 = plt.subplot(212, sharex=ax1)
plt.imshow(ordered_locus, interpolation='nearest', cmap='coolwarm')
plt.show()

def show_correlation(idx):
    re_fltr=np.logical_not(np.logical_or(np.isnan(flds[idx,:]),
                                       np.any(np.isnan(ordered_locus), axis=0)))

    pre_fields = flds[idx, re_fltr]
    pf = pre_fields.reshape(1, pre_fields.shape[0])
    sp = ordered_locus[:, re_fltr]

    ax1 = plt.subplot(211)
    plt.imshow(pf,
               interpolation='nearest', cmap='coolwarm')
    plt.setp(ax1.get_xticklabels(), fontsize=6)
    ax2 = plt.subplot(212, sharex=ax1)
    plt.imshow(sp,
               interpolation='nearest', cmap='coolwarm')
    plt.show()

# show_correlation(1)
rho, p_val = spearmanr(ordered_locus, axis=1)

fltr = np.logical_not(np.logical_and(p_val < 0.05, np.absolute(rho) > 0.3))
p_val[fltr] = np.NaN
rho[fltr] = np.NaN
np.fill_diagonal(p_val, np.NaN)
np.fill_diagonal(rho, np.NaN)


plt.subplot(211)
plt.imshow(rho,
           interpolation='nearest', cmap='coolwarm')
plt.colorbar()
plt.subplot(212)
plt.imshow(p_val,
           interpolation='nearest', cmap='coolwarm')
plt.colorbar()
plt.show()

def chr_correlation(idex):
    re_fltr=np.logical_not(np.logical_or(np.isnan(flds[idex,:]),
                                   np.any(np.isnan(ordered_locus), axis=0)))
    pre_fields = flds[idex, re_fltr]
    pf = pre_fields.reshape(1, pre_fields.shape[0])
    sp = ordered_locus[:, re_fltr]
    np.apply_along_axis(spearmanr, 1, sp, pf.flatten())
    intermediate = np.apply_along_axis(spearmanr, 1, sp, pf.flatten())
    return np.split(intermediate.T, 2)

def gross_accumulate():
    rho = []
    p_val = []
    for i in range(0, flds.shape[0]):
        t_rho, t_p_val = chr_correlation(i)
        rho.append(t_rho.flatten())
        p_val.append(t_p_val.flatten())
    return np.array(rho), np.array(p_val)

t_rho, t_p_val = gross_accumulate()

selector = np.logical_not(np.logical_and(t_p_val < 0.05, np.absolute(t_rho) > 0.3))

t_rho[selector] = np.NaN
t_p_val[selector] =np.NaN

drugnames = np.array([drug for index, drug in sorted(drug_idx_rv.items())])

sel2 = np.any(np.logical_not(np.isnan(t_rho)), axis=1)
t_rho = t_rho[sel2,:]
t_p_val = t_p_val[sel2,:]
drugnames = drugnames[sel2].tolist()

show_matrix_with_names(t_rho,
                       drugnames,
                       range(1, 25),
                       colormap = 'coolwarm')

show_matrix_with_names(t_p_val,
                       drugnames,
                       range(1, 25),
                       colormap = 'coolwarm')

show_correlation(drug_idx['Epirubicin'])

show_correlation(drug_idx['Bortezomib'])
