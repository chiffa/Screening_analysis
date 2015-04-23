__author__ = 'ank@stowers.org'

import os
import xlrd
import numpy as np
from collections import defaultdict
from csv import reader, writer
from pprint import pprint
from time import sleep
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, minimize
from itertools import combinations

debug = False
aspects = {
        'debug':['flatten_and_group.pre_process']
        }

source_folder = 'L:\\ank\\real-runs\\'

past_mappings = os.path.join(source_folder, 'conditions-to-include-in-clustering_2010-09-16.csv') # used by Pavelka in his clustering
relative_growth = os.path.join(source_folder, 'meanNormalizedFinalAverageSpotIntensity.csv')
spot_locations = os.path.join(source_folder, 'spotLocation.csv')
output_location = os.path.join(source_folder, 'Re_analysis_ank-2.csv')
total_log = os.path.join(source_folder, 'Re_analysis_ank-failed_to_parse.csv')
output_folder = os.path.join(source_folder, 'Re_analysis_ank-2')


def growth(timepoints, stepness, maxVal, midpoint, delay):
    # maxVal = 50  # works well
    # we will just remember all the parameters; this functions seems to be working better for some reason then the other ones.
    if delay > 0: # delay control
        timepoints = timepoints - delay
        timepoints[timepoints < 0] = 0
    return maxVal/(1. + np.exp(-np.log(2)*stepness*(timepoints - midpoint)))  # TODO: problem, lag and midpoint are redundant


def minimize_function_builder(timepoints, means, errors):

    def minfunct(paramset):
        stepness, maxVal, midpoint, delay = paramset
        estimate = growth(timepoints, stepness, maxVal, midpoint, delay)
        difference = np.abs(estimate - means)
        # ponderated_difference = difference/np.abs(errors)
        ponderated_difference = difference
        return np.sum(ponderated_difference**2)

    return minfunct


def rec_path_join(base, *args):
    for elt in args:
        base = os.path.join(base, elt)
    return base


def build_unified_mappings(mappings_file):
    storage_dict = defaultdict(list)
    authorised_replicates = {}
    with open(mappings_file, 'rb') as source:
        csv_reader = reader(source)
        csv_reader.next()
        for line in csv_reader:
            if line:
                standard_name = line[0]
                standard_folder = rec_path_join(source_folder, line[1], line[2])
                storage_dict[standard_name].append(standard_folder)
                authorised_replicates[standard_folder] = line[3].split(';')

    storage_dict = dict([(key, tuple(value)) for key, value in storage_dict.iteritems()])
    return storage_dict, authorised_replicates


def check_unified_mappings(unified_mapping, authorised_replicas):
    for name, location_list in unified_mapping.iteritems():
        for location in location_list:
            if not os.path.isdir(location):
                print 'location %s for %s is mapped wrong' % (location, name)


# read location-depended spot locations:
def build_spot_map():
    spot_dict = defaultdict(dict)
    with open(spot_locations, 'rb') as source:
        csv_reader = reader(source)
        for line in csv_reader:
            if line[3]:
                spot_dict[line[1]].update({line[2]: line[3]})
    return spot_dict


def read_file(file_path, spots_map):
    with open(file_path, 'rb') as source:
        csv_reader = reader(source)
        names =[]
        times = np.array(csv_reader.next()[1:])
        table = []
        for line in csv_reader:
            if line[0] in spots_map.keys():
                names.append(spots_map[line[0]])
                table.append(line[1:])

        return np.array(names), times.astype(np.float), np.array(table).astype(np.float)


# walks the source folder and reads the original data.
def pull_curves(name2folder, folder2replicas, spots_map):
    name2full_location = defaultdict(list)
    for name, folders in name2folder.iteritems():
        for folder in folders:
            for subfolder in os.walk(folder):
                parsed_subfolder = subfolder[0].split('\\')
                if parsed_subfolder[-1][-2:] in ['_A', '_B', '_C'] \
                        and os.path.isdir(subfolder[0]):
                    for fle in subfolder[2]:
                        if 'averageSpotIntensity' in fle and parsed_subfolder[-1][-1] in folder2replicas[folder]:
                            name2full_location[name].append((parsed_subfolder[-1][-1], os.path.join(subfolder[0], fle)))

    read_out_map = defaultdict(list)
    # dico[condition] = [replica_1, replica_2, ....]
    # where replica_x = (aneuploid indexes, time_index, data_matrix)
    for name, fle_list in name2full_location.iteritems():
        for replicate, fle_to_read in fle_list:
            current_map = spots_map[replicate]
            current_dataframe = read_file(fle_to_read, current_map)
            read_out_map[name].append(current_dataframe)

    return read_out_map


def flatten_and_group(condition_specific_replica, condition_name):

    def show(subdict, name):
        plt.title(name)
        for time, value_list in subdict.iteritems():
            for value in value_list:
                plt.plot(time, value)
        plt.show()


    def high_level_pre_process():
        # a totally different solution: fit each replicate individually,
        # output individual fit repeats; average parameters => In addition to failure to raise, also gives a lag value
        pass


    def pre_process(time, value_list):
        # set every lane to 0
        value_list = np.array(value_list)
        value_list = value_list - np.min(np.mean(value_list, axis=0))

        take_off = np.max(np.min(value_list, axis=0)[:-1])
        if take_off < 5:
            return  value_list, 'fail on take-off: %s; ' % take_off, 1

        mns = np.mean(value_list, axis=0)
        if np.any(mns > 10):
            abs_diff = np.max(value_list, axis=0) - np.min(value_list, axis=0)
            rel_diffs = abs_diff / mns
            divergence = np.max(abs_diff[mns > 10])
            relative_divegence = np.max(rel_diffs[mns > 10])
            if divergence > 15 or relative_divegence > 0.9:
                fail = 'fail on repeats: %s ; %s; '%(divergence, relative_divegence)
                return  value_list, fail, 2

        return value_list, '', 0


    def process(time, value_list):
        # Combine for each element all the x and y columns by concatenating them; then plot
        # Use regression to fit in a curve
        newdict = {}
        for value in value_list:
            plt.plot(time, value)
        value_list = np.array(value_list)
        mns = np.mean(value_list, axis=0)
        errs = np.std(value_list, axis=0, ddof=1)
        try:
            popt, pcov = curve_fit(growth, np.array(time), mns, p0=[0.1, 60, 60, 40], sigma=errs)
            # print 'popt:', popt
            # print 'pcov:', pcov
            fit = (growth(np.array(time), *popt),
                   [1./popt[0]] + popt[1:].tolist(),
                   np.mean(np.abs(growth(np.array(time), *popt)-mns)))
        except RuntimeError:
            # print 'fit not found'
            fit = (mns, [-1., -1., -1.], -1.)
        newdict[time] = (mns, errs, fit) #fit-related piping. Maybe a Monad would be better here?
        print fit[1], fit[2]
        return newdict, []


    def process2(time, value_list):
        # Combine for each element all the x and y columns by concatenating them; then plot
        # Use regression to fit in a curve
        newdict = {}
        for time, value_list in subdict.iteritems():
            for value in value_list:
                plt.plot(time, value)
            value_list = np.array(value_list)
            mns = np.mean(value_list, axis=0)
            errs = np.std(value_list, axis=0, ddof=1)
            mfunct = minimize_function_builder(np.array(time), mns, errs)
            OR_object = minimize(mfunct, [0.1, 60, 60, 40], bounds=[(0.01, 0.9), (40, 100), (0, 500), (0, 150)])
            popt = OR_object.x
            # print OR_object.x
            # print OR_object.message
            # print 'expected error:', mfunct(popt)
            if OR_object.success:
                fit = (growth(np.array(time), *popt),
                       [1./popt[0]] + popt[1:].tolist(),
                       np.mean(np.abs(growth(np.array(time), *popt)-mns)))
            else:
                # print OR_object.message
                fit = (mns,
                       [-1., -1., -1., -1.],
                       -1.)
            newdict[time] = (mns, errs, fit) #fit-related piping. Maybe a Monad would be better here?
        return newdict, fit[1] + [fit[2]]

    def show_processed(subdict, name):
        plt.title(name)
        for time, (mns, stds, fit) in subdict.iteritems():
            plt.errorbar(time, mns, stds, fmt='k.')
            plt.plot(time, fit[0], 'k', label=' doubling time: %.2f h\n max: %.2f \n midpoint: %.0f h\n lag: %.0f h\n error: %.2f\n '% tuple(fit[1] + [fit[2]]))
            plt.legend(loc='upper left', prop={'size':10})
        plt.show()


    aneuploid2MergedTable = defaultdict(lambda : [[], []])

    for replicate in condition_specific_replica:
        for i, aneuploid_ID in enumerate(replicate[0]):
            time = tuple(replicate[1].tolist())
            value = replicate[2][i, :].tolist()
            aneuploid2MergedTable[aneuploid_ID][0] = time
            aneuploid2MergedTable[aneuploid_ID][1].append(value)

    supercollector = []
    fail_collector = []
    for aneuploid_ID, (time, value_list) in aneuploid2MergedTable.iteritems():
        # run the pre-procssing round
        re_value_list, fail_message, errvalue = pre_process(time, value_list)
        # record all the raw data!!!!!
        for repeat in value_list:
            fail_collector.append([aneuploid_ID, condition_name, errvalue] + repeat)
        # if it just fails to rise or the fit is unperfect, record it anyway
        if errvalue < 2:
            if errvalue == 0:
                re_dict, carryover = process2(time, re_value_list)
            else:
                re_dict, carryover = ({}, ['inf', 'NA', 'NA', 'NA', np.mean(np.abs(np.mean(re_value_list, axis=0)))])
            supercollector.append([aneuploid_ID, condition_name] + carryover)

    return supercollector, fail_collector


def iterate_through_conditions(readout_map):
    super_collector = [['strain', 'condition', 'stepness', 'maxVal', 'midpoint', 'delay',  'fit error']]
    fail_collector = []
    for condition, condition_specific_replica in readout_map.iteritems():
        d_super, d_fail = flatten_and_group(condition_specific_replica, condition)
        super_collector += d_super
        fail_collector += d_fail

    return super_collector, fail_collector


# finally, write out the resulting curves to a destination file
def write_out_curves(locations, out_path):
    with open(out_path, 'wb') as source:
        csv_writer = writer(source)
        for line in locations:
            csv_writer.writerow(line)


if __name__ == "__main__":
    canonical_mappings, canonical_replicas = build_unified_mappings(past_mappings)
    check_unified_mappings(canonical_mappings, canonical_replicas)
    # pprint(dict(canonical_mappings))
    spots_dict = build_spot_map()
    readout_map = pull_curves(canonical_mappings, canonical_replicas, spots_dict)
    collector, fails = iterate_through_conditions(readout_map)
    print np.array(collector)
    write_out_curves(collector, output_location)
    write_out_curves(fails, total_log)
