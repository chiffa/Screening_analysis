from matplotlib import pyplot as plt
from csv import reader as csv_reader
import numpy as np
import os
import matplotlib
from chiffatools.linalg_routines import gini_coeff, rm_nans

source_folder = 'L:\\Users\\andrei\\Jin_2010_re_analysis'
values = 'speeds_v.csv'
errors = 'speeds_err.csv'


def csv2numpy(source, c_header=True, r_header=True):

    def correct_line(_row):
        return [float(item) if item not in ['inf', '', ' '] else np.inf for item in _row]

    with open(source, 'r') as source_file:
        reader = csv_reader(source_file)

        if c_header:
            c_headers = reader.next()
        else:
            c_headers = []

        r_headers = []
        data_container = []

        for row in reader:
            if r_header:
                r_headers.append(row[0])
                row = row[1:]
            data_container.append(correct_line(row))

        return np.array(data_container), c_headers, r_headers


def plot_err_bar_plot():
    vals, c_headers, r_headers = csv2numpy(os.path.join(source_folder, values))
    errs, _, _ = csv2numpy(os.path.join(source_folder, errors))

    c_headers = c_headers[1:]
    r_headers = r_headers
    print vals.shape
    print errs.shape
    print len(c_headers)
    print len(r_headers)
    ginis = []
    averages = []

    x_ladder = np.arange(len(c_headers))
    for i in range(0, len(r_headers)):

        rev_vals = 1./vals[i, :]
        support = np.sum(rev_vals != 0)
        # rev_vals = rev_vals[rev_vals != 0]

        gini = gini_coeff(rev_vals)
        ginis.append(gini)
        mean = np.mean(rev_vals)
        averages.append(mean)

        # vals_only = vals[i, :]
        # vals_only = vals_only[vals_only != np.inf]
        # gini = gini_coeff(vals_only)
        # mean = np.mean(vals_only)

        plt.errorbar(x_ladder, vals[i, :],
                     yerr=errs[i, :],
                     label='%s - gini %.2f - mean %.2f - support %s' % (r_headers[i],
                                                                        gini, mean,
                                                                        support),
                     marker='o',
                     ls='None')

    plt.xticks(x_ladder, c_headers, rotation='vertical')
    plt.legend()
    # matplotlib.rcParams.update({'font.size': 8})
    plt.show()

    r_headers = np.array(r_headers)
    ginis = np.array(ginis)
    plt.plot(ginis, averages, 'ko', label='aneuploids')
    averages = np.array(averages)

    print r_headers

    index = [42, 46]
    plt.plot(ginis[index], averages[index], 'ro', label=', '.join(r_headers[index]))

    index = [43, 45]
    plt.plot(ginis[index], averages[index], 'go', label=', '.join(r_headers[index]))

    index = [44, 47]
    plt.plot(ginis[index], averages[index], 'bo', label=', '.join(r_headers[index]))

    plt.legend()
    plt.ylabel('average fitness')
    plt.xlabel('gini index')
    plt.show()

    index = 42
    ladder = np.sort(1./vals[index, :])
    plt.plot(x_ladder, np.cumsum(ladder)/np.sum(ladder), 'b', label='generalist - gini %.2f' %
                                                                    gini_coeff(ladder))
    plt.plot(x_ladder, np.linspace(0, 1., len(c_headers)), 'r')

    index = 45
    ladder = np.sort(1./vals[index, :])
    plt.plot(x_ladder, np.cumsum(ladder)/np.sum(ladder), 'k', label='specialist - gini %.2f' %
                                                                    gini_coeff(ladder))

    index = 37
    ladder = np.sort(1./vals[index, :])
    plt.plot(x_ladder, np.cumsum(ladder)/np.sum(ladder), 'm', label='strong specialist - gini '
                                                                    '%.2f' % gini_coeff(ladder))

    plt.legend()
    plt.ylabel('cumulative fitness share')
    plt.xlabel('stresses')

    plt.show()


if __name__ == "__main__":
    plot_err_bar_plot()
