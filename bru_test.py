__author__ = 'ank'

from scipy.sparse import lil_matrix
import numpy as np
from scipy.linalg import cho_factor, cho_solve#

fudge = 10e-8


def build_IO_currents_array(IO_index_pair, shape):
    IO_array = np.zeros((shape[0], 1))
    IO_array[IO_index_pair[0], 0], IO_array[IO_index_pair[1], 0] = (1.0, -1.0)
    return  IO_array


def get_voltages(conductivity_laplacian, IO_index_pair):
    # dense_laplacian = conductivity_laplacian + fudge*np.identity(conductivity_laplacian.shape[0])
    # dense_laplacian = conductivity_laplacian + np.ones(conductivity_laplacian.shape)
    # v1 = np.ones((conductivity_laplacian.shape[0], 1))
    v1 = np.array([[1, 0, 0, 0]]).T
    v2 = np.array([[1, 0, 0, 0]])
    dense_laplacian = conductivity_laplacian + np.dot(v1, v2)
    print dense_laplacian
    IO_array = build_IO_currents_array(IO_index_pair, conductivity_laplacian.shape)
    chdec = cho_factor(dense_laplacian)
    voltages = cho_solve(chdec,IO_array)
    return '{0:.4f}'.format(float(abs(voltages[IO_index_pair[0]]-voltages[IO_index_pair[1]])*3/5)**2*5)



conductivity_laplacian = np.array([ [2, -1, -1, 0],
                                    [-1, 2, 0, -1],
                                    [-1, 0, 2, -1],
                                    [0, -1, -1, 2]  ])

print conductivity_laplacian

print get_voltages(conductivity_laplacian, (0, 3))