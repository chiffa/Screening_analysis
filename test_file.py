__author__ = 'ank'

import numpy as np

nparr = np.array(['1','2','3.14','1e-3','b','nan','inf','-inf'])

nparr = nparr.reshape((2, 4))

print type(nparr), nparr.shape
print nparr
print np.genfromtxt(nparr)