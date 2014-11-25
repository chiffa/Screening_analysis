__author__ = 'ank'

import numpy as np
from matplotlib import pyplot as plt
from itertools import combinations
from chiffatools. Linalg_routines import rm_nans
from scipy.stats import ttest_ind


def t_test_matrix(current_lane, breakpoints):
    dst = len(breakpoints)
    p_vals = np.empty((dst, dst))
    p_vals.fill(np.NaN)
    subsets = np.split(current_lane, breakpoints[:-1])
    for i, j in combinations(range(0, dst), 2):
        _, p_val = ttest_ind(rm_nans(subsets[i]), rm_nans(subsets[j]), False)
        p_vals[i, j] = p_val
    return p_vals


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def pull_breakpoints(contingency_list):
    no_nans_parsed = rm_nans(contingency_list)
    contingency = np.lib.pad(no_nans_parsed[:-1] == no_nans_parsed[1:], (1, 0), 'constant', constant_values=(True, True))
    nans_contingency = np.zeros(contingency_list.shape).astype(np.bool)
    nans_contingency[np.logical_not(np.isnan(contingency_list))] = contingency
    breakpoints = np.nonzero(np.logical_not(nans_contingency))[0].tolist()
    return breakpoints


def show_breakpoints(breakpoints):
    for point in breakpoints:
        plt.axvline(x=point, color='b')


def generate_breakpoint_mask(breakpoints):
    support = np.zeros((np.max(breakpoints), ))
    pre_brp = 0
    for i, brp in enumerate(breakpoints):
        support[pre_brp:brp] = i
        pre_brp = brp
    return support


def inflate_support(length, breakpoints, values=None):
    if values is None:
        values = np.array(range(0, len(breakpoints)))
    ret_array = np.zeros((50, length))
    for _i in range(1, values.shape[0]):
        ret_array[:, breakpoints[_i-1]: breakpoints[_i]] = values[_i]
    return ret_array