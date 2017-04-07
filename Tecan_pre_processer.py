import os
import xlrd
import numpy as np
from itertools import product
from Support import debug_wrapper
from matplotlib import pyplot as plt
import matplotlib as mlb
from scipy.ndimage.filters import gaussian_filter1d as gf_1d
from collections import defaultdict
from csv import reader, writer
from pprint import pprint
from time import time
from string import ascii_lowercase
from sklearn.gaussian_process import GaussianProcess
from multiprocessing import Pool, current_process
from chiffatools.dataviz import smooth_histogram

tinit = time()
mlb.rcParams['font.size'] = 10.0
mlb.rcParams['figure.figsize'] = (26, 15)

# todo: add a map of the wells and perform a statistical significance analysis
# todo: add discounting for the yeast division lag in the new cells.

# file_location = 'U:/ank/2015/TcanScreen/03.26.2015/GertonLabTcan/'
file_location = 'C:\\Users\\Andrei\\Desktop'
# file_name = 'Tecan_9-26-2014.xlsx'
# file_name = 'Book1.xlsx'
file_name = '11232016growth.xlsx'
d_time = 30./60.

time_map = defaultdict(int)
current_path = os.path.join(file_location, file_name)
wb = xlrd.open_workbook(current_path)


@debug_wrapper
def extract_plate(a1_coordinates, sheet):
    plate = np.zeros((8, 12))
    for i, j in product(range(0, 8), range(0, 12)):
        _i, _j = (np.array([i, j]) + np.array(list(a1_coordinates))).astype(np.uint32).tolist()
        # print _i, _j, sheet.cell(_i, _j).value
        if sheet.cell(_i, _j).value:
            plate[i, j] = sheet.cell(_i, _j).value
        else:
            print 'missing data @', _i, _j
            plate[i, j] = -1
    return plate


def extract_plate_dit():
    collect_dict = []
    for s in wb.sheets():
        if 'Magellan' in s.name:
            latest = 0
            for row in range(1, s.nrows):
                if row % 9 == 1:
                    plate = extract_plate((row, 1), s)
                    collect_dict.append(plate)
                    latest += d_time
    return np.array(collect_dict)


def smooth_plate(plate, window):
    re_plate = plate.copy()
    for i, j in product(range(0, 8), range(0, 12)):
        re_plate[:, i, j] = gf_1d(plate[:, i, j], window)
    return re_plate


def plot_growth(plates_stack, grad=False):
    fig = plt.gcf()
    fig.canvas.set_window_title('min:%s, max:%s, total points: %s'%(np.min(plates_stack), np.max(plates_stack), plates_stack.shape[0]))
    for i, j in product(range(0, 8), range(0, 12)):
        fig = plt.subplot(8, 12, i*12+j+1)
        data = plates_stack[:, i, j]
        if grad:
            plt.title('%s, %s'%('{0:.0f}'.format(float(d_time/np.max(data)*60)), '{0:.2f}'.format(d_time*np.argmax(data)-(d_time/np.max(data)*3))))
        plt.plot(data.tolist())
        plt.ylim((np.min(plates_stack), np.max(plates_stack)))
        if j != 0:
            fig.set_yticklabels([])
        if i != 7:
            fig.set_xticklabels([])
        if i == 7:
            tick_lbls = ['{0:.1f}'.format(d_time*int(item)) for item in fig.get_xticks()]
            # tick_lbls = [(lambda x, term: term if x % 2 == 0 else '')(_i, val) for _i, val in enumerate(tick_lbls)]
            fig.set_xticklabels(tick_lbls)
            for tick in fig.xaxis.get_major_ticks():
                tick.label.set_fontsize(10)
                tick.label.set_rotation('vertical')

    plt.savefig('im_%s.png'%int(time()-tinit), dpi=300)
    plt.show()
    # plt.clf()


def group_plot(plates_stack, zoomlist):
    timepad = np.linspace(0, d_time*plates_stack.shape[0], num=plates_stack.shape[0])
    for sublist in zoomlist:
        legend = []
        for elt in sublist:
            plt.plot(timepad, plates_stack[:,elt[0], elt[1]])
            legend.append(str('%s %s')%(ascii_lowercase[elt[0]], elt[1]+1))
        plt.legend(legend, loc='upper left')
        plt.savefig('group_im_%s.png'%int(time()-tinit))
        plt.show()
        # plt.clf()


def analyse(plates_stack, zoomlist):
    if intermediate_show:
        plot_growth(plates_stack)
    reference_std = np.std(plates_stack[:, 0, 0])*2
    print reference_std
    log_stack = np.log10(plates_stack)/np.log10(2)
    for i, j in product(range(0, 8), range(0, 12)):
        log_stack[:, i, j] = log_stack[:, i, j] - np.mean(log_stack[range(0, 3), i, j])
    if intermediate_show:
        plot_growth(log_stack)
    grad_stack = np.zeros(log_stack.shape)
    for i, j in product(range(0, 8), range(0, 12)):
        grad_stack[:, i, j] = np.gradient(gf_1d(log_stack[:, i, j], 2))
    if intermediate_show:
        plot_growth(grad_stack, True)
    # group_plot(plates_stack, zoomlist)
    group_plot(log_stack, zoomlist)
    # group_plot(grad_stack, zoomlist)


def correct(plate, position, injections):
    new_plate = np.zeros((plate.shape[0]+injections-1, plate.shape[1], plate.shape[2]))
    new_plate[:position+1, :, :] = plate[:position+1, :, :]
    new_plate[position+injections-1:, :, :] = plate[position:, :, :]
    diffplate = (plate[position+1, :, :] - plate[position, :, :]) / float(injections)
    for i in range(1, injections):
        new_plate[position+i, :, :] = plate[position, :, :] + diffplate * i
    return new_plate


def gaussian_process_regress(timeseries, std, timestamps=None, show=False):

    def show_routine():
        fig = plt.figure()
        plt.errorbar(timestamps.ravel(), timeseries, errors, fmt='r.', markersize=10, label=u'Observations')
        plt.plot(pre_timestamps, y_pred, 'b-', label=u'Prediction')
        plt.fill(np.concatenate([pre_timestamps, pre_timestamps[::-1]]),
                np.concatenate([y_pred - 1.9600 * sigma,
                               (y_pred + 1.9600 * sigma)[::-1]]),
                alpha=.5, fc='b', ec='None', label='95% confidence interval')
        plt.xlabel('$x$')
        plt.ylabel('$f(x)$')
        plt.legend(loc='upper left')

        plt.show()

    if timestamps is None:
        timestamps = np.linspace(0, timeseries.shape[0]*d_time, timeseries.shape[0])[:, np.newaxis]

    pre_timestamps = timestamps.copy()

    keep_mask = timeseries > 0.0001
    timestamps = timestamps[:, 0][keep_mask][:, np.newaxis]
    timeseries = timeseries[keep_mask]

    nugget = np.convolve(timeseries, np.ones((5,))/5, mode='valid')
    nugget = np.lib.pad(nugget, (2, 2), 'edge')
    errors = np.abs(nugget - timeseries)
    errors[errors < std] = std
    nugget = np.power((errors), 2)

    gp = GaussianProcess(regr='linear', corr='squared_exponential', theta0=10,
                     thetaL=1e-1, thetaU=100,
                     nugget=nugget,
                     random_start=100)

    gp.fit(timestamps, timeseries)
    y_pred, MSE = gp.predict(pre_timestamps, eval_MSE=True)
    sigma = np.sqrt(MSE)

    if show:
        show_routine()

    elif np.any(y_pred < 0.00001):
        # show_routine()
        pass

    return y_pred, sigma


def gaussian_process_wrapper(bulk_arguments):
    i, j, pl, std = bulk_arguments
    print 'loess', i, j
    return ((i,j), gaussian_process_regress(pl, std))


def map_adapter(plate, std):
     for i, j in product(range(0, 8), range(0, 12)):
         yield i,j, plate[:, i, j], std


def loess(plate):
    refsample = plate_3D_array[generate_reference_mask(plate)]
    std = np.std(refsample)

    plate = plate - np.percentile(refsample[refsample > 0.0001], 0.5)
    re_plate = plate.copy()
    fine_tune = np.percentile(refsample[refsample > 0.0001], 1)

    retset = map(gaussian_process_wrapper, map_adapter(plate, std))
    for ((i, j), (ret, _)) in retset:
        re_plate[:, i, j] = ret
    re_plate[re_plate < fine_tune] = fine_tune
    return re_plate


def generate_reference_mask(plate):

    def extracted_growth_detector(_1D_array):
        return (np.percentile(_1D_array, 97.5) - np.percentile(_1D_array, 2.5)) < 0.10

    ref_mask = np.zeros((8, 12)).astype(np.bool)
    ref_mask[:, 0] = True
    ref_mask[:, 11] = True
    ref_mask[0, :] = True
    ref_mask[7, :] = True

    growth_detected = np.apply_along_axis(extracted_growth_detector, 0, plate)
    ref_mask = np.logical_and(ref_mask, growth_detected)
    timed_ref_mask = np.repeat(ref_mask[np.newaxis, :, :], plate.shape[0], axis=0)

    return timed_ref_mask


def smooth_and_interpolate(plate):
    refpoints = loess(plate)
    # deviation = np.abs(refpoints - plate)
    # _99 = np.percentile(deviation, 99)
    # plate[deviation > _99] = -1
    # refpoints = loess(plate)
    return refpoints


def del_exception(plate, position):
    new_plate = np.zeros((plate.shape[0]-1, plate.shape[1], plate.shape[2]))
    new_plate[:position, :, :] = plate[:position, :, :]
    print position
    new_plate[position:, :, :] = plate[position+1:, :, :]
    return new_plate


def del_range(plate, positionList):
    for position in positionList:
        plate = del_exception(plate, position)
    return plate


def L2_dist(plates_stack, growth_times):
    pass

def L2_growth_difference(growth_speed_matrix):
    pass

# TODO: write repeat similarity distance => L2 on growth curves; L2 on growth conditions
# TODO: return the test as two matrices
# TODO: perform the matrix hierarchical clustering


if __name__ == "__main__":
    intermediate_show = True
    plate_3D_array = extract_plate_dit()
    plate_3D_array = smooth_and_interpolate(plate_3D_array)

    # plate_3D_array = smooth_plate(plate_3D_array, 2)
    # plate_3D_array = del_exception(plate_3D_array, 220)
    # plate_3D_array = del_exception(plate_3D_array, 220)
    # plate_3D_array = correct(plate_3D_array, 219, 6)
    # del_range(plate_3D_array, range(220,222))
    zoomlist = []
    zoomlist = [
                # [(1, 1), (5, 1), (1, 3), (4, 3), (5, 3) ],
                # [(1, 1), (1, 3), (1, 6), (1, 10)],
                # [(3, 1), (3, 3), (3, 6), (3, 10)],
                # [(5, 1), (4, 3), (5, 3)],
                # [(2, 6), (6, 6), (3, 7), (4, 7), (1, 2), (2, 2)],
                ]
    analyse(plate_3D_array, zoomlist)
