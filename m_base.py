__author__ = 'ank'

from itertools import combinations
from matplotlib import pyplot as plt

Bp_Dict=[
        230218,
        813184,
        316620,
        1531933,
        576874,
        270161,
        1090940,
        562643,
        439888,
        745751,
        666816,
        1078177,
        924431,
        784333,
        1091291,
        948066,
]

print Bp_Dict, len(Bp_Dict)

rg = range(0,len(Bp_Dict))

appPlod_list = []
AllList = []
tl = float(sum(Bp_Dict))

for i in range(0,len(rg)):
    appPlod_list = appPlod_list + [ len(chrlist)/16.0+1 for chrlist in combinations(rg, i)]

for i in range(0,len(rg)):
    for chrset in combinations(rg,i):
        AllList.append(sum(Bp_Dict[i] for i in chrset)/tl+1)


bins = [ i/16.0+1 for i in range(0, 16, 1)]
plt.hist(appPlod_list, bins=bins, facecolor='b', alpha=0.5, histtype='step', normed=True)
# plt.show()
plt.hist(AllList, bins=bins, facecolor='r', alpha=0.5, histtype='step', normed=True)
# plt.normpdf( bins, 8, 2)
plt.show()

# plt.hist([appPlod_list, AllList], bins=16, histtype='bar', color=['crimson', 'chartreuse'],
#                             label=['Real Ploidy', 'Apparent Ploidy'], normed=True)
# plt.legend()
# plt.show()
