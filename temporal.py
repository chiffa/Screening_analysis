__author__ = 'ank'

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from pylab import get_cmap

def show_matrix_with_names(matrix, vert_names, horiz_names, colormap='b_jet', overlay=None):

    def full_rename(namelist, subtypes):

        def rename(list, subtype, character):
            ret_list = []
            for bearer in list:
                modded = False
                for mod_decider in subtype:
                    if mod_decider in bearer:
                        modded = True
                        print bearer
                        print character
                        ret_list.append(bearer+' '+str(character))
                        break
                if not modded:
                    ret_list.append(bearer)
            return ret_list

        new_vert_names = vert_names
        for idux, subtype in enumerate(vert_names):
            if idux > 0:
                print ''.join(['*']*idux)
                new_vert_names = rename(new_vert_names, subtype, ''.join(['*']*idux))

        return new_vert_names



    if colormap == 'b_jet':
        prism_cmap = get_cmap('jet')
        prism_vals = prism_cmap(np.arange(0, 1, 0.01))
        prism_vals[99] = [0, 0, 0, 1]
        costum_cmap = colors.LinearSegmentedColormap.from_list('my_colormap', prism_vals)

    else:
        costum_cmap = get_cmap(colormap)

    plt.imshow(matrix, interpolation='nearest', cmap=costum_cmap)
    plt.colorbar()
    if overlay:
        # we need to transform out of imshow
        plt.scatter(overlay[0][0], overlay[0][1], label=overlay[2])
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand", borderaxespad=0.)
    plt.tick_params(axis='both', labelsize=10)

    if type(vert_names) == list:
        plt.yticks(range(0, len(vert_names)), vert_names, rotation='horizontal')
    if type(horiz_names) == list:
        plt.xticks(range(0, len(horiz_names)), horiz_names, rotation='vertical')
    if type(vert_names) == tuple:
        vert_names = full_rename(vert_names[0], vert_names[1:])
        plt.yticks(range(0, len(vert_names)), vert_names, rotation='horizontal')
    if type(horiz_names) == tuple:
        horiz_names = full_rename(horiz_names[0], horiz_names[1:])
        plt.xticks(range(0, len(horiz_names)), horiz_names, rotation='vertical')


    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.show()