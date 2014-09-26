__author__ = 'ank'

import os
import xlrd
import numpy as np
from itertools import product
from Support import debug_wrapper
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d as gf_1d
from collections import defaultdict
from csv import reader, writer
from pprint import pprint

file_location = 'U:/ank/2014/Screen_with_Jin'
file_name = 'Tecan_9-26-2014.xlsx'
d_time = 15./60.

current_path = os.path.join(file_location, file_name)
wb = xlrd.open_workbook(current_path)


@debug_wrapper
def extract_plate(a1_coordinates, sheet):
    plate = np.zeros((8, 12))
    for i, j in product(range(0, 8), range(0, 12)):
        _i, _j = (np.array([i, j]) + np.array(list(a1_coordinates))).astype(np.uint32).tolist()
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
    plt.figure(figsize=(30.0, 20.0))
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
    if grad:
        plt.savefig('grads.png', dpi=500)
    plt.show()
    plt.clf()

def analyse(plates_stack):
    plot_growth(plates_stack)
    reference_std = np.std(plates_stack[:, 0, 0])*2
    print reference_std
    log_stack = np.log10(plates_stack)/np.log10(2)
    plot_growth(log_stack)
    grad_stack = np.zeros(log_stack.shape)
    for i, j in product(range(0, 8), range(0, 12)):
        grad_stack[:, i, j] = np.gradient(gf_1d(log_stack[:, i, j], 1.5))
    plot_growth(grad_stack, True)
    #TODO: add discounting for the yeast difvision lag in the new cells.

if __name__ == "__main__":
    plate_3D_array = extract_plate_dit()
    analyse(plate_3D_array)