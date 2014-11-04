__author__ = 'ank'

from csv import reader, writer
import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Put the name of the source file containing as the first row the names of genes and as the following rows the growth rate
# of 1N, 2N and 3N colonies. The file has to be a comma-separated file (.csv). The extension doesn't matter
loc = 'C:\\Users\\ank\\Desktop\\hj-clustering.csv'
# Where the method would write out the assigned clusters
loc2 = 'C:\\Users\\ank\\Desktop\\hj-clustering-out.csv'

# these are the parameters that define the clustering behavior. Please modify them at will.
epsilon = 0.5
minimal_samples = 20


# Following lines import the data and format it as required for modules
namearray = []
tuple_array = []
trituple_array = []


with open(loc,'rb') as source:
    rd = reader(source)
    rd.next()
    for row in rd:
        namearray.append(row[0])
        tuple_array.append([float(elt) for elt in row[1:-1]])
        trituple_array.append([float(elt) for elt in row[1:]])


namearray = np.array(namearray)
tuple_array = np.array(tuple_array)
trituple_array = np.array(trituple_array)
x,y,z = trituple_array.T.tolist()

print namearray, tuple_array, trituple_array

##################################################################################
# Uncomment the section below for a 2-D clustering
##################################################################################

db = DBSCAN(eps=epsilon, min_samples=minimal_samples).fit(tuple_array)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
lbls = db.labels_

print lbls

unique_labels = set(lbls)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))

print type(lbls)


for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (lbls == k)

    xy = tuple_array[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)

    xy = tuple_array[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6)

plt.show()

##################################################################################
# Uncomment the section below for a 3-D clustering
##################################################################################

# db = DBSCAN(eps=epsilon, min_samples=minimal_samples).fit(trituple_array)
# core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
# core_samples_mask[db.core_sample_indices_] = True
# lbls = db.labels_
# unique_labels = set(lbls)
# colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# for k, col in zip(unique_labels, colors):
#     if k == -1:
#         # Black used for noise.
#         col = 'k'
#
#     class_member_mask = (lbls == k)
#
#     xy = trituple_array[class_member_mask & core_samples_mask]
#     ax.scatter(xy[:, 0], xy[:, 1], xy[:, 2], c=col, s=30)
#
#     xy = trituple_array[class_member_mask & ~core_samples_mask]
#     ax.scatter(xy[:, 0], xy[:, 1], xy[:, 2], c=col, s=20)
#
# ax.set_xlabel('1 N')
# ax.set_ylabel('2 N')
# ax.set_zlabel('3 N')
#
# plt.show()

#############################################################################################################
# This is the method that will be used for result output indepently of the 2 or 3-D clustering
#############################################################################################################

farr = np.array([namearray, trituple_array.T[0], trituple_array.T[1], trituple_array.T[2] ,lbls]).T.tolist()
print farr

with open(loc2,'wb') as source:
    wr = writer(source)
    for row in farr:
        wr.writerow(row)