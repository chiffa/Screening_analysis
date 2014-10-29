__author__ = 'ank'

import numpy as np
from matplotlib import pyplot as plt

debug = 0

def rs(matrix, name):
    plt.title(name)
    plt.imshow(matrix, interpolation='nearest')
    plt.colorbar()
    plt.show()
    plt.clf()

def debug_wrapper(funct):

    def check_matrix(*args,**kwargs):
        result = funct(*args, **kwargs)
        if debug:
            rs(result, funct.__name__)
        return result

    return check_matrix

def pretty_gradual_plot(data, concentrations, strain_name_map, drug_name, blank_line=200):

    def inner_scatter_plot(mean, std, relative, limiter=4):
        series = np.zeros(mean.shape)
        cell_type = np.zeros(mean.shape)
        for i, name in enumerate(names):
            series[i, :] = np.arange(i, c.shape[0]*(len(names)+40)+i, len(names)+40)
            cell_type[i, :] = i
            plt.scatter(series[i, :], mean[i, :], c=cm(i/float(len(names))), s=35, label=name)
        plt.errorbar(series.flatten(), mean.flatten(), yerr=std.flatten(), fmt=None, capsize=0)
        plt.xticks(np.mean(series, axis=0), c)
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=len(names)/limiter, mode="expand", borderaxespad=0.,prop={'size':6})
        if not relative:
            plt.axhline(y=blank_line)
        plt.show()

    filter = np.all(np.logical_not(np.isnan(data)), axis=(1, 2))
    names = [strain_name_map[i] for i in filter.nonzero()[0].tolist()]
    c = concentrations[filter, :][0, :]
    mean = np.nanmean(data[filter, :, :], axis=-1)
    std = np.nanstd(data[filter, :, :], axis=-1)
    cm = plt.cm.get_cmap('spectral')

    refmean = mean[:, 0].reshape((mean.shape[0], 1))
    refstd = std[:, 0].reshape((mean.shape[0], 1))
    rel_mean, rel_std = (mean/refmean, np.sqrt(np.power(refstd, 2)+np.power(std, 2))/mean)

    inner_scatter_plot(mean, std, False)
    inner_scatter_plot(rel_mean, rel_std, True)

    mean_mean = np.nanmean(mean, axis=0)
    std_mean = np.nanstd(mean, axis=0)
    mean_std = np.nanmean(std, axis=0)
    total_std = np.sqrt(np.power(std_mean, 2) + np.power(mean_std, 2))
    confusables = np.sum(mean - std < blank_line, axis=0) / float(len(names))

    rel_mean_mean = np.nanmean(rel_mean, axis=0)
    rel_std_mean = np.nanstd(rel_mean, axis=0)
    rel_mean_std = np.nanmean(rel_std, axis=0)
    rel_total_std = np.sqrt(np.power(rel_std_mean, 2) + np.power(rel_mean_std, 2))

    plt.subplot(212)
    plt.plot(mean_mean, c=cm(0), label='mean of mean')
    plt.plot(mean_std, c=cm(.25), label='mean of std')
    plt.plot(std_mean, c=cm(.50), label='std of mean')
    plt.plot(total_std, c=cm(0.75), label='total std')
    plt.axhline(y=blank_line)

    plt.subplot(211)
    plt.plot(rel_mean_mean, c=cm(0), label='mean of mean')
    plt.plot(rel_mean_std, c=cm(.25), label='mean of std')
    plt.plot(rel_std_mean, c=cm(.50), label='std of mean')
    plt.plot(rel_total_std, c=cm(0.75), label='total std')
    plt.plot(confusables, c=cm(0.90), label='confusable with null')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand", borderaxespad=0.,prop={'size':8})

    plt.show()


def pretty_get_resistant_susceptible(data, strain_name_map, drug_name, blank_line=200):

    def inner_scatter_plot(mean, std, names, selector=None, relative=False, limiter=4):
        if selector is not None:
            mean = mean[selector, :]
            std = std[selector, :]
            names = np.array(names)[selector].tolist()
        series = np.zeros(mean.shape)
        cell_type = np.zeros(mean.shape)
        for i, name in enumerate(names):
            series[i, :] = np.arange(i, mean.shape[1]*(len(names)+40)+i, len(names)+40)
            cell_type[i, :] = i
            plt.scatter(series[i, :], mean[i, :], c=cm(i/float(len(names))), s=35, label=name)
        plt.errorbar(series.flatten(), mean.flatten(), yerr=std.flatten(), fmt=None, capsize=0)
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=len(names)/limiter, mode="expand", borderaxespad=0.,prop={'size':6})
        if not relative:
            plt.axhline(y=blank_line)
        plt.show()

    filter = np.all(np.logical_not(np.isnan(data)), axis=(1, 2))
    reverse_filter = filter.nonzero()[0]
    names = [strain_name_map[i] for i in filter.nonzero()[0].tolist()]
    mean = np.nanmean(data[filter, :, :], axis=-1)
    std = np.nanstd(data[filter, :, :], axis=-1)
    cm = plt.cm.get_cmap('spectral')

    refmean = mean[:, 0].reshape((mean.shape[0], 1))
    refstd = std[:, 0].reshape((mean.shape[0], 1))
    rel_mean, rel_std = (mean/refmean, np.sqrt(np.power(refstd, 2)+np.power(std, 2))/mean)

    # From here, we are processing just the best survivors v.s. worst survivors
    dead = mean-std<blank_line

    # first, filter out all those that cannot be distinguished from background
    upper = np.percentile(rel_mean, 95, axis = 0).reshape(1, rel_mean.shape[1])
    lower = np.percentile(rel_mean, 5, axis = 0).reshape(1, rel_mean.shape[1])
    upper[0, :-2] = 100
    lower[0, :-2] = 0

    destroyed = np.any(dead, axis=1).nonzero()[0]
    resistant = list(set(np.any(rel_mean>upper, axis=1).nonzero()[0])-set(destroyed))
    susceptible = list(set(np.any(rel_mean<lower, axis=1).nonzero()[0])-set(destroyed))

    ret_array = np.empty(data.shape[0])
    ret_array.fill(np.NaN) # we might try to do somethign to distinguish the "no effect" from "no data"
    ret_array[reverse_filter[resistant]] = 1
    ret_array[reverse_filter[susceptible]] = 0
    ret_array[reverse_filter[destroyed]] = -1

    if len(resistant) < 1:
        int_ret = pretty_get_resistant_susceptible(data[:,:-1,:], strain_name_map, drug_name, blank_line=blank_line)
        int_ret[int_ret>0] = int_ret[int_ret>0] * 0.5
        flt = int_ret > 0
        ret_array[flt] = int_ret[flt]
        return ret_array

    inner_scatter_plot(mean, std, names, resistant, False, 1)
    inner_scatter_plot(rel_mean, rel_std, names, resistant, True, 1)

    inner_scatter_plot(mean, std, names, susceptible, False, 1)
    inner_scatter_plot(rel_mean, rel_std, names, susceptible, True, 1)

    inner_scatter_plot(mean, std, names, destroyed, False, 1)
    inner_scatter_plot(rel_mean, rel_std, names, destroyed, True, 1)

    print 'resistant:', np.array(names)[resistant].tolist()
    print 'susceptible:', np.array(names)[susceptible].tolist()
    print 'destroyed:', np.array(names)[destroyed].tolist()

    return ret_array

def get_resistant_susceptible(data, drug_name, blank_line=200):

    filter = np.all(np.logical_not(np.isnan(data)), axis=(1, 2))
    reverse_filter = filter.nonzero()[0]
    mean = np.nanmean(data[filter, :, :], axis=-1)
    std = np.nanstd(data[filter, :, :], axis=-1)

    refmean = mean[:, 0].reshape((mean.shape[0], 1))
    refstd = std[:, 0].reshape((mean.shape[0], 1))
    rel_mean, rel_std = (mean/refmean, np.sqrt(np.power(refstd, 2)+np.power(std, 2))/mean)

    # From here, we are processing just the best survivors v.s. worst survivors
    dead = mean-std<blank_line

    # first, filter out all those that cannot be distinguished from background
    upper = np.percentile(rel_mean, 95, axis = 0).reshape(1, rel_mean.shape[1])
    lower = np.percentile(rel_mean, 5, axis = 0).reshape(1, rel_mean.shape[1])
    upper[0, :-2] = 100
    lower[0, :-2] = 0

    destroyed = np.any(dead, axis=1).nonzero()[0]
    resistant = list(set(np.any(rel_mean>upper, axis=1).nonzero()[0])-set(destroyed))
    susceptible = list(set(np.any(rel_mean<lower, axis=1).nonzero()[0])-set(destroyed))

    ret_array = np.empty(data.shape[0])
    ret_array.fill(np.NaN)
    ret_array[reverse_filter] = 0
    ret_array[reverse_filter[resistant]] = 1
    ret_array[reverse_filter[susceptible]] = -0.5
    ret_array[reverse_filter[destroyed]] = -1

    if len(resistant)<1:
        int_ret = get_resistant_susceptible(data[:,:-1,:], drug_name, blank_line=blank_line)
        int_ret[int_ret>0] = int_ret[int_ret>0]*0.5
        flt = int_ret > 0
        ret_array[flt] = int_ret[flt]

    return ret_array