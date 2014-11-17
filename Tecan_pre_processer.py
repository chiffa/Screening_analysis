__author__ = 'ank'

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

tinit = time()
mlb.rcParams['font.size'] = 10.0
mlb.rcParams['figure.figsize'] = (30,20)
# todo: add discounting for the yeast division lag in the new cells.
# todo: Normalize OD with respect to the empty wells

file_location = 'U:/ank/2014/Screen_with_Jin/11.04.2014'
# file_name = 'Tecan_9-26-2014.xlsx'
file_name = 'Book1.xlsx'
d_time = 15./60.


current_path = os.path.join(file_location, file_name)
wb = xlrd.open_workbook(current_path)


@debug_wrapper
def extract_plate(a1_coordinates, sheet):
    plate = np.zeros((8, 12))
    for i, j in product(range(0, 8), range(0, 12)):
        _i, _j = (np.array([i, j]) + np.array(list(a1_coordinates))).astype(np.uint32).tolist()
        # print _i, _j, sheet.cell(_i, _j).value
        plate[i, j] = sheet.cell(_i, _j).value
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


def plot_growth(plates_stack, grad=False):
    fig = plt.gcf()
    fig.canvas.set_window_title('min:%s, max:%s, total points: %s'%(np.min(plates_stack), np.max(plates_stack), plates_stack.shape[0]))
    for i, j in product(range(0, 8), range(0, 12)):
        fig = plt.subplot(8, 12, i*12+j+1)
        data = plates_stack[:, i, j]
        if grad:
            plt.title('%s, %s'%('{0:.2f}'.format(float(d_time/np.max(data))), d_time*int(np.argmax(data))))
        plt.plot(data.tolist())
        plt.ylim((np.min(plates_stack), np.max(plates_stack)))
        fig.set_xticklabels([])
        fig.set_yticklabels([])
    plt.savefig('im_%s.png'%int(time()-tinit), dpi=500)
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
    plot_growth(plates_stack)
    reference_std = np.std(plates_stack[:, 0, 0])*2
    print reference_std
    log_stack = np.log10(plates_stack)/np.log10(2)
    for i,j in product(range(0, 8), range(0, 12)):
        log_stack[:, i, j] = log_stack[:, i, j] - np.mean(log_stack[range(0, 3), i, j])
    plot_growth(log_stack)
    grad_stack = np.zeros(log_stack.shape)
    for i, j in product(range(0, 8), range(0, 12)):
        grad_stack[:, i, j] = np.gradient(gf_1d(log_stack[:, i, j], 2))
    plot_growth(grad_stack, True)
    group_plot(plates_stack, zoomlist)
    # group_plot(log_stack, zoomlist)
    group_plot(grad_stack, zoomlist)


def correct(plate, position, injections):
    new_plate = np.zeros((plate.shape[0]+injections-1, plate.shape[1], plate.shape[2]))
    new_plate[:position+1, :, :] = plate[:position+1, :, :]
    new_plate[position+injections-1:, :, :] = plate[position:, :, :]
    diffplate = (plate[position+1, :, :] - plate[position, :, :]) / float(injections)
    for i in range(1, injections):
        new_plate[position+i, :, :] = plate[position, :, :] + diffplate * i
    return new_plate


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


if __name__ == "__main__":
    plate_3D_array = extract_plate_dit()
    # plate_3D_array = del_exception(plate_3D_array, 220)
    # plate_3D_array = del_exception(plate_3D_array, 220)
    # plate_3D_array = correct(plate_3D_array, 219, 6)
    # del_range(plate_3D_array, range(220,222))
    plate_3D_array = plate_3D_array - np.min(plate_3D_array) + 0.001
    zoomlist = []
    # zoomlist = [
    #             [(5, 6), (6, 6), (5, 2), (6, 2)],
    #             [(5, 6), (6, 6), (3, 6), (4, 6)],
    #             [(5, 6), (6, 6), (1, 6), (2, 6)],
    #             [(5, 6), (6, 6), (3, 7), (4, 7)],
    #             [(2, 6), (6, 6), (3, 7), (4, 7), (1, 2), (2, 2)],
    #             ]
    analyse(plate_3D_array, zoomlist)
