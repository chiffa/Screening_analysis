"""
In it's current form performs an analysis of assay results by fitting a saturated exponential
growth curve to it
"""
import os
import numpy as np
from collections import defaultdict
from csv import reader, writer
from matplotlib import pyplot as plt
from scipy.optimize import minimize as sc_min
from pickle import dump, load
import pandas as pd
import seaborn.apionly as sns
from chiffatools.linalg_routines import rm_nans, gini_coeff
import matplotlib as mlb
from pprint import pprint

mlb.rcParams['font.size'] = 10.0
mlb.rcParams['figure.figsize'] = (25, 15)


debug = False
aspects = {
        'debug': ['flatten_and_group.pre_process']
        }

# source_folder = 'L:\\Users\\andrei\\real-runs\\'
source_folder = 'C:\\Users\\Andrei\\Documents\\real-runs'

past_mappings = os.path.join(source_folder, 'conditions-to-include-in-clustering_2010-09-16.csv')
# conditions and replicates used by Pavelka in his clustering
past_mappings = os.path.join(source_folder, 'conditions-to-include-in-clustering_reduced.csv') # Nature paperclustering
relative_growth = os.path.join(source_folder, 'meanNormalizedFinalAverageSpotIntensity.csv')
spot_locations = os.path.join(source_folder, 'spotLocation.csv')
output_location = os.path.join(source_folder, 'Re_analysis_ank-3.csv')
total_log = os.path.join(source_folder, 'Re_analysis_ank-failed_to_parse.csv')
output_folder = os.path.join(source_folder, 'Re_analysis_ank-3')


def growth(time_points, steepness, max_val, midpoint, offset):
    # offset can be used as pre-drop
    preliminary = max_val / (1. + np.exp(-np.log(2) * steepness * (time_points - midpoint))) + offset
    return preliminary


def build_fit_quality_estimator(time_points, od_means):

    def fit_quality_estimator(parameter_set):
        steepness, max_val, midpoint, delay = parameter_set
        estimate = growth(time_points, steepness, max_val, midpoint, delay)
        difference = np.abs(estimate - od_means)
        return np.sum(difference**2)

    return fit_quality_estimator


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
        names = []
        times = np.array(csv_reader.next()[1:])
        table = []
        for line in csv_reader:
            if line[0] in spots_map.keys():
                names.append(spots_map[line[0]])
                table.append(line[1:])

        return np.array(names), times.astype(np.float), np.array(table).astype(np.float)


# walks the source folder and reads the original data.
def pull_curves(name2folder, folder_2_replicas, spots_map):

    def errf(x):
        raise x


    name2full_location = defaultdict(list)
    for name, folders in name2folder.iteritems():
        for folder in folders:
            for subfolder in os.walk(folder, onerror=errf):
                parsed_subfolder = subfolder[0].split('\\')
                if parsed_subfolder[-1][-2:] in ['_A', '_B', '_C'] \
                        and os.path.isdir(subfolder[0]):
                    for fle in subfolder[2]:
                        if 'averageSpotIntensity' in fle and parsed_subfolder[-1][-1] in folder_2_replicas[folder]:
                            print 'appending location', os.path.join(subfolder[0], fle)
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


def merge_replicates(condition_specific_replicas):
    aneuploid_2_merged_replicates = defaultdict(lambda: [[], []])

    for replicate in condition_specific_replicas:
        for i, aneuploid_id in enumerate(replicate[0]):
            time = tuple(replicate[1].tolist())
            value = replicate[2][i, :].tolist()

            if not aneuploid_2_merged_replicates[aneuploid_id][0]:
                aneuploid_2_merged_replicates[aneuploid_id][0] = time

            if time == aneuploid_2_merged_replicates[aneuploid_id][0]:
                aneuploid_2_merged_replicates[aneuploid_id][1].append(value)

    return aneuploid_2_merged_replicates


def flatten_and_group(condition_specific_replica, condition_name, aneuploid_index):

    def nested_list_from_len(def_integer):
        return [[] for _ in range(0, def_integer)]

    def nan_3diff(val1, val2):
        if val1 in ['NA', 'inf'] or val2 in ['NA', 'inf']:
            return 'NA'
        else:
            return val1 - 6.6*val2

    def nan_3sum(val1, val2):
        if val1 in ['NA', 'inf'] or val2 in ['NA', 'inf']:
            return 'NA'
        else:
            return val1 + 6.6*val2

    def iterative_fit(time_points, od_means):

        def random_between(min_val, max_val):
            return np.random.uniform(min_val, max_val)

        od_mean_arr = np.array(od_means)
        od_mean_arr -= np.min(od_mean_arr)

        bounds = [(0.01, 0.25), (10, 100), (0, 200), (0, 1)]

        parameters, errcode = fit_parameters(time_points, od_mean_arr, params_bounds=bounds)
        if parameters[-1] > 1 and errcode != 1:
            for i in range(0, 5):
                start = [random_between(*bound) for bound in bounds]
                parameters, errcode = fit_parameters(time_points, od_mean_arr,
                                                     starting_params=start, params_bounds=bounds)
                if parameters[-1] <= 1:
                    break

            if parameters[-1] > 1.5:  # error too big
                errcode = 2

        return od_mean_arr, parameters, errcode

    def fit_parameters(time_points, od_means, starting_params=[0.16, 50., 60., 0.],
                       params_bounds=[(0.05, 0.5), (10, 150), (0, 200), (-2, 2)]):

        take_off = np.max(od_means[1:-1])

        if take_off < 10:
            return ['inf', 'NA', 'NA', 'NA', np.mean(np.abs(np.mean(od_means, axis=0)))], 1

        fit_quality_estimator = build_fit_quality_estimator(np.array(time_points), od_means)
        optimizer = sc_min(fit_quality_estimator, starting_params,
                           method='L-BFGS-B', bounds=params_bounds)
        optimal_parameters = optimizer.x

        if optimizer.success:
            rise_start = nan_3diff(optimal_parameters[2], 1./optimal_parameters[0])
            rise_end = nan_3sum(optimal_parameters[2], 1./optimal_parameters[0])
            err_no = 0

            if rise_end > time_points[-1]:
                # print 'lag too great - timeline runs out before fit is reached'
                err_no = 3

            time_points = np.array(time_points)
            mask = np.logical_and(time_points < rise_end, time_points > rise_start)

            if np.sum(mask.astype(np.int8)) < 5:
                print 'too few points in the growth area to perform a fit'
                err_no = 4

            return [1./optimal_parameters[0]] + optimal_parameters[1:].tolist() \
                        + [np.mean(np.abs(growth(np.array(time_points),
                                                 *optimal_parameters) - od_means))],\
                        err_no

        else:
            print optimizer.message,
            print optimizer.x
            optimal_parameters = optimizer.x
            return [1./optimal_parameters[0]] + optimal_parameters[1:].tolist() +\
                        [np.mean(np.abs(growth(np.array(time_points),
                                               *optimal_parameters) - od_means))],\
                        -1

    def show(time_points, value, fit_parameters, name):
        plt.title(name)
        time_points = np.array(time_points)
        higher_time = np.linspace(np.min(time_points), np.max(time_points), 100)
        plt.plot(time_points, value, 'ro')
        plt.plot(higher_time, growth(higher_time, 1 / fit_parameters[0], *fit_parameters[1:-1]),
                 'k', label=' doubling time_points: %.2f h\n max: %.2f \n midpoint: %.0f h\n '
                            'offset: %.2f \n error: %.2f\n exp-phase: %.2f h - %.2f h\n ' % tuple(
                fit_parameters+[nan_3diff(fit_parameters[2], fit_parameters[0]),
                                nan_3sum(fit_parameters[2], fit_parameters[0])]))
        plt.legend(loc='upper left', prop={'size': 10})
        plt.show()

    fit_param_collector = []
    od_means_collector = []
    base_shape = len(aneuploid_index.keys())
    doubling_time_lane = nested_list_from_len(base_shape)
    midpoints_lane = nested_list_from_len(base_shape)
    delay_lane = nested_list_from_len(base_shape)
    maxval_lane = nested_list_from_len(base_shape)
    lastval_lane = nested_list_from_len(base_shape)

    aneuploid_2_merged_replicates = merge_replicates(condition_specific_replica)

    for aneuploid_ID, (time_points, od_means_list) in aneuploid_2_merged_replicates.iteritems():
        print aneuploid_ID,
        doubling_time_holder = []
        midpoint_holder = []
        delay_lane_holder = []
        maxval_lane_holder = []
        lastval_lane_holder = []
        for od_mean in od_means_list:
            norm_repeat, fit_params, error_code = iterative_fit(time_points, od_mean)

            if error_code in [2, 4]:

                fit_params = [np.nan, np.nan, np.nan, np.nan, np.mean(np.abs(np.mean(od_mean, axis=0)))]

            # if error_code in [0, 2, 3, 4] and '16C(2)' in condition_name:
            #     show(time_points, norm_repeat, fit_params, aneuploid_ID+', '+condition_name)
            # if error_code == 2 or error_code == -1:
            #     show(time_points, norm_repeat, fit_params, aneuploid_ID+', '+condition_name)

            od_means_collector.append([aneuploid_ID, condition_name, error_code] + od_mean)
            fit_param_collector.append([aneuploid_ID, condition_name, error_code] + fit_params)
            a_i = aneuploid_index[aneuploid_ID]
            doubling_time_holder.append(fit_params[0])
            midpoint_holder.append(fit_params[2])
            delay_lane_holder.append(nan_3diff(fit_params[2], fit_params[0]))
            maxval_lane_holder.append(fit_params[1])
            lastval_lane_holder.append(od_mean[-1])

        doubling_time_lane[a_i] = doubling_time_holder
        midpoints_lane[a_i] = midpoint_holder
        delay_lane[a_i] = delay_lane_holder
        maxval_lane[a_i] = maxval_lane_holder
        print lastval_lane_holder
        lastval_lane[a_i] = lastval_lane_holder

    return fit_param_collector, od_means_collector, doubling_time_lane, midpoints_lane, \
           delay_lane, maxval_lane, lastval_lane


def iterate_through_conditions(readout_map):
    fit_params_collector = [['strain', 'condition', 'fitting result', 'doubling time(h)', 'maxVal',
                        'midpoint', 'delay',  'fit error']]
    od_means_collector = []
    condition_names = []
    speeds = []
    midpoints = []
    delays = []
    maxvals = []
    lastvals = []

    aneuploids_set = set()

    for _, condition_specific_replicas in readout_map.iteritems():
        for replica in condition_specific_replicas:
            aneuploids_set.update(set(replica[0]))

    aneuploids_names = sorted(list(aneuploids_set))
    aneuploids_index = dict([(aneup_name, _i) for _i, aneup_name in enumerate(aneuploids_names)])

    for condition, condition_specific_replicas in readout_map.iteritems():
        print '\n\n', condition
        condition_names += [condition]
        fit_params_coll, od_means_coll, doubling_time_lane, midpoints_lane, delay_lane, \
        maxval_lane, lastval_lane = \
            flatten_and_group(condition_specific_replicas, condition, aneuploids_index)
        fit_params_collector += fit_params_coll
        od_means_collector += od_means_coll
        speeds.append(doubling_time_lane)
        midpoints.append(midpoints_lane)
        delays.append(delay_lane)
        maxvals.append(maxval_lane)
        lastvals.append(lastval_lane)
        print '\n'
        pprint(lastval_lane)

    cons_obj = (aneuploids_names, condition_names, speeds, midpoints, delays, maxvals, lastvals)
    return fit_params_collector, od_means_collector, cons_obj


# finally, write out the resulting curves to a destination file
def write_out_curves(matrix, out_path):
    with open(out_path, 'wb') as source:
        csv_writer = writer(source)
        for line in matrix:
            csv_writer.writerow(line)


def reduce_table(conservation_object):
    # TODO: looks that it's here that the lastvalues fail.
    def reduction_routine(list_to_reduce):
        list_to_reduce = [str(elt) for elt in list_to_reduce]
        num_list = np.genfromtxt(np.array(list_to_reduce))
        non_numerical = np.logical_or(np.isnan(num_list), np.logical_not(np.isfinite(num_list)))
        if sum(non_numerical) == len(list_to_reduce):
            return np.inf, np.nan
        if sum(non_numerical) == len(list_to_reduce)-1:
            return np.nan, np.nan
        numerical_redux = num_list[np.logical_not(non_numerical)]
        mn, sd = (np.nanmean(numerical_redux),
                  1.96*np.nanstd(numerical_redux, ddof=1)/np.sqrt(len(list_to_reduce)-sum(
                      non_numerical)))
        if sd/mn < 0.5 and mn > 2:  # noise acceptable for a large enough mean
            return mn, sd
        if mn < 2:  # noise does not matter if we are close to 0
            return mn, sd

        else:
            return np.nan, np.nan

    def high_order_reduction(embedded_list):
        return np.array([[reduction_routine(lst) for lst in cond_lst] for cond_lst in embedded_list])

    def split_pandas_frame(data):
        # print data.shape
        df_v = pd.DataFrame(data[:, :, 0].T, aneup_names, condition_names)
        df_err = pd.DataFrame(data[:, :, 1].T, aneup_names, condition_names)
        return df_v, df_err

    def render(variable1, variable2, name):
        plt.subplot(1, 2, 1)
        plt.title(name+' values')
        sns.heatmap(variable1)
        plt.subplot(1, 2, 2)
        plt.title(name+' errors')
        sns.heatmap(variable2)
        plt.show()

    def errplot_with_selectors(table, errtable):
        for i in range(0, len(selector)):
            v1 = table.reset_index().values[:, 1:][i, :].flatten()
            v2 = errtable.reset_index().values[:, 1:][i, :].flatten()
            v1 = v1.astype(np.float)
            nl = np.sum(np.logical_not(np.isnan(v1)))
            gni = gini_coeff(1./rm_nans(v1))
            plt.errorbar(condition_index, v1, v2, fmt='.', label='%s; gini: %.2f, valid: %s' %
                                                                 (pre_selector[i], gni, nl))
        plt.xticks(condition_index, condition_names, rotation='vertical')
        plt.gca().set_yscale("log", nonposy='clip')
        plt.legend(loc='upper left', prop={'size': 10})
        plt.show()

    def ratio_with_errs(v1, err1, v2, err2):
        verification_array = np.array((v1, v2, err1, err2))
        if np.all(np.logical_and(np.logical_not(np.isnan(verification_array)),
                                 np.isfinite(verification_array))):
            return v2 / v1, np.sqrt((v2/v1)**2*(err1**2/v1**2+err2**2/v2**2))
            # ATTENTION: since we are interested in
            # growth speed ratios, we are using an inverted ratio. Thus v_ref/v_calc is perfectly normal
        else:
            return np.nan, np.nan

    def ratio_wrapper(v1_vect, err1_vect, v2_vect, err2_vect):
        return [ratio_with_errs(v1_vect[_i], err1_vect[_i], v2_vect[_i], err2_vect[_i])
                for _i in range(0, len(v1_vect))]

    def diff_with_errs(v1, err1, v2, err2):
        verification_array = np.array((v1, v2, err1, err2))
        if np.all(np.logical_and(np.logical_not(np.isnan(verification_array)), np.isfinite(verification_array))):
            return v1 - v2, np.sqrt(err1**2 + err2**2)
        else:
            return np.nan, np.nan

    # numpy vectorize?
    def diff_wrapper(v1_vect, err1_vect, v2_vect, err2_vect):
        return [diff_with_errs(v1_vect[_i], err1_vect[_i], v2_vect[_i], err2_vect[_i])
                for _i in range(0, len(v1_vect))]

    def spin_conditions():
        for condition in condition_names:
            lag = np.array(lag_v.loc[:, condition]).astype(np.float).flatten()
            lag_e = np.array(lag_errs.loc[:, condition]).astype(np.float).flatten()
            ratio = np.array(twister_v.loc[:, condition]).astype(np.float).flatten()
            ratio_e = np.array(twister_errs.loc[:, condition]).astype(np.float).flatten()
            plt.title(condition)
            for _j, a_name in enumerate(aneup_names):
                if np.isnan(lag[_j]) or not np.isfinite(lag[_j]):
                    plt.errorbar(lag[_j], ratio[_j], ratio_e[_j], lag_e[_j], fmt='.w', label=a_name)
                else:
                    plt.annotate(a_name, (lag[_j], ratio[_j]))
                    if np.abs(lag[_j]-1) < 5 and np.abs(ratio[_j]-1) < 0.1:
                        plt.errorbar(lag[_j], ratio[_j], ratio_e[_j], lag_e[_j],
                                     fmt='.k', label=a_name)
                    elif lag[_j]>10 and ratio[_j] > 1:
                        plt.errorbar(lag[_j], ratio[_j], ratio_e[_j], lag_e[_j],
                                     fmt='.r', label=a_name)
                    else:
                        plt.errorbar(lag[_j], ratio[_j], ratio_e[_j], lag_e[_j],
                                     fmt='.b', label=a_name)
            plt.legend(loc='upper left', prop={'size': 10})
            plt.savefig('%s.png' % condition)
            # plt.show()
            plt.clf()

    aneup_names, condition_names, speeds, midpoints, delays, maxvals, lastvals = conservation_object

    speeds_v, speeds_err = split_pandas_frame(high_order_reduction(speeds))
    midpoints_v, midpoints_err = split_pandas_frame(high_order_reduction(midpoints))
    delays_v, delays_err = split_pandas_frame(high_order_reduction(delays))
    maxvals_v, maxvals_err = split_pandas_frame(high_order_reduction(maxvals))
    lastvals_v, lastvals_err = split_pandas_frame(high_order_reduction(lastvals))

    pprint(lastvals)

    pprint(lastvals_v)

    raw_input('press enter to continue')

    speeds_v.to_csv('speeds_v.csv')
    speeds_err.to_csv('speeds_err.csv')
    midpoints_v.to_csv('midpoints_v.csv')
    midpoints_err.to_csv('midpoints_err.csv')
    delays_v.to_csv('delays_v.csv')
    delays_err.to_csv('delays_err.csv')
    maxvals_v.to_csv('maxvals_v.csv')
    maxvals_err.to_csv('maxvals_err.csv')
    lastvals_v.to_csv('lastvals_v.csv')
    lastvals_err.to_csv('lastvals_err.csv')

    re_index = dict([(name, _i) for _i, name in enumerate(aneup_names)])

    pre_selector = ['U1', 'U2', 'U3']
    pre_selector = aneup_names

    selector = np.array([re_index[pre_s] for pre_s in pre_selector])
    condition_index = np.array([i for i, _ in enumerate(condition_names)])
    # bigger_index = np.repeat(condition_index, len(selector))
    # v1 = speeds_v.reset_index().values[:, 1:][selector, :].flatten()
    # v2 = speeds_err.reset_index().values[:, 1:][selector, :].flatten()

    # errplot_with_selectors(speeds_v, speeds_err)
    # errplot_with_selectors(delays_v, delays_err)

    s_v = speeds_v.reset_index().values[:, 1:]
    s_err = speeds_err.reset_index().values[:, 1:]

    ref_v = s_v[re_index['U1'], :]
    ref_err = s_err[re_index['U1'], :]

    twister_fused = np.array([ratio_wrapper(s_vs, s_errs, ref_v, ref_err)
                              for s_vs, s_errs in zip(s_v.tolist(), s_err.tolist())])
    twister_fused = np.rollaxis(twister_fused, 1)
    twister_v, twister_errs = split_pandas_frame(twister_fused)

    l_v = delays_v.reset_index().values[:, 1:]
    l_err = delays_err.reset_index().values[:, 1:]

    ref_v = l_v[re_index['U1'], :]
    ref_err = l_err[re_index['U1'], :]

    lag_fused = np.array([diff_wrapper(l_vs, l_errs, ref_v, ref_err)
                          for l_vs, l_errs in zip(l_v.tolist(), l_err.tolist())])
    lag_fused = np.rollaxis(lag_fused, 1)
    lag_v, lag_errs = split_pandas_frame(lag_fused)

    # print lag_v
    # print lag_errs

    # TODO: recalculate the euploid from the aneuploids

    spin_conditions()

    # errplot_with_selectors(twister_v, twister_errs)

    # find the one we actually want by a 2_d plot

    # ix_1 = re_index['controlHaploid']
    # ix_2 = re_index['U1']
    #
    # print np.array(twister_errs.iloc[[ix_2]]).astype(np.float)
    #
    # plt.errorbar(np.array(twister_v.iloc[[ix_1]]).astype(np.float).flatten(),
    #              np.array(twister_v.iloc[[ix_2]]).astype(np.float).flatten(),
    #              xerr=np.array(twister_errs.iloc[[ix_1]]).astype(np.float).flatten(),
    #              yerr=np.array(twister_errs.iloc[[ix_2]]).astype(np.float).flatten())
    # plt.show()

    twister_v.to_csv('twister_v.csv')
    twister_errs.to_csv('twister_errs.csv')
    #
    (twister_v - twister_errs).to_csv('twister_v_conservative.csv')

    # render(speeds_v, speeds_err, 'duplication_time')
    # render(midpoints_v, midpoints_err, 'midpoints')
    # render(delays_v, delays_err, 'delays')
    # render(np.log10(twister_v), twister_errs, 'relative_growth_speed')
    # render(lag_v, lag_errs, 'relative lags')


def regress_euploid(conservation_object):

    def reduction_routine(list_to_reduce):
        list_to_reduce = [str(elt) for elt in list_to_reduce]
        num_list = np.genfromtxt(np.array(list_to_reduce))
        non_numerical = np.logical_or(np.isnan(num_list), np.logical_not(np.isfinite(num_list)))
        if sum(non_numerical) == len(list_to_reduce):
            return 0, np.nan
        if sum(non_numerical) == len(list_to_reduce)-1:
            return np.nan, np.nan
        numerical_redux = 1. / num_list[np.logical_not(non_numerical)]
        mn, sd = (np.nanmean(numerical_redux), 1.96*np.nanstd(numerical_redux, ddof=1)/np.sqrt(len(
            list_to_reduce)-sum(non_numerical)))
        if np.log10(mn) > np.log10(sd) and mn < 0.5:
            return mn, sd
        else:
            return np.nan, np.nan

    def errplot_with_selectors(table, errtable):
        for i in range(0, len(selector)):
            v1 = table.reset_index().values[:, 1:][selector[i], :].flatten()
            v2 = errtable.reset_index().values[:, 1:][selector[i], :].flatten()
            v1 = v1.astype(np.float)
            nl = np.sum(np.logical_not(np.isnan(v1)))
            gni = gini_coeff(v1)
            plt.errorbar(condition_index, v1, v2, fmt='.', label='%s; gini: %.2f, valid: %s' %
                                                                 (pre_selector[i], gni, nl))
        plt.xticks(condition_index, condition_names, rotation='vertical')
        plt.gca().set_yscale("log", nonposy='clip')
        plt.legend(loc='upper left', prop={'size':10})
        plt.show()

    def higher_reduction(embedded_list):
        return np.array([[reduction_routine(lst) for lst in cond_lst] for cond_lst in embedded_list])

    def split_pandas_frame(data):
        df_v = pd.DataFrame(data[:, :, 0].T, aneup_names, condition_names)
        df_err = pd.DataFrame(data[:, :, 1].T, aneup_names, condition_names)
        return df_v, df_err

    aneup_names, condition_names, speeds, midpoints, delays = conservation_object
    speeds_v, speeds_err = split_pandas_frame(higher_reduction(speeds))
    midpoints_v, midpoints_err = split_pandas_frame(higher_reduction(midpoints))
    delays_v, delays_err = split_pandas_frame(higher_reduction(delays))

    re_index = dict([(name, _i) for _i, name in enumerate(aneup_names)])

    # pre_selector = ['U1', 'U2', 'U3']
    pre_selector = aneup_names[:-7]
    selector = np.array([re_index[pre_s] for pre_s in pre_selector])
    condition_index = np.array([(i) for i, _ in enumerate(condition_names)])

    current_table = speeds_v.reset_index().values[:, 1:][selector, :].astype(np.float)
    current_table_err = speeds_err.reset_index().values[:, 1:][selector, :].astype(np.float)
    aneuploid_gini_indexes = np.apply_along_axis(gini_coeff, 1, current_table)
    aneuploid_mean_survival = np.apply_along_axis(np.nanmean, 1, current_table)

    # print aneuploid_mean_survival

    argsort_indexes = np.argsort(aneuploid_gini_indexes)

    # print aneuploid_gini_indexes[argsort_indexes]

    pre_selector = ['U1']
    selector = np.array([re_index[pre_s] for pre_s in pre_selector])
    l_idx = selector[-1]

    # print np.nanmax(current_table)

    for j, i in enumerate(range(4, 7)):
        euploid_reconstruction = np.apply_along_axis(np.nanmean, 0, current_table[argsort_indexes[:i], :])
        euploid_reconstruction = euploid_reconstruction / np.nanmean(euploid_reconstruction) * np.nanmax(aneuploid_mean_survival)+0.001
        euploid_err_reconstruction = np.apply_along_axis(np.nanmean, 0, current_table_err[argsort_indexes[:i], :])

        # print i, gini_coeff(euploid_reconstruction), np.nanmean(euploid_reconstruction),\
            # np.nansum((euploid_reconstruction - speeds_err.reset_index().values[:, 1:][l_idx, :].astype(np.float))**2)

        speeds_v.loc['reconstruction %s' % i] = euploid_reconstruction
        speeds_err.loc['reconstruction %s' % i] = euploid_err_reconstruction
        pre_selector.append('reconstruction %s' % i)
        aneup_names.append('reconstruction %s' % i)

    # pre_selector = aneup_names
    re_index = dict([(name, _i) for _i, name in enumerate(aneup_names)])
    selector = np.array([re_index[pre_s] for pre_s in pre_selector])
    selector = np.array([re_index[pre_s] for pre_s in pre_selector])
    # errplot_with_selectors(speeds_v, speeds_err)


if __name__ == "__main__":
    canonical_mappings, canonical_replicas = build_unified_mappings(past_mappings)
    check_unified_mappings(canonical_mappings, canonical_replicas)
    # pprint(dict(canonical_mappings))
    spots_dict = build_spot_map()
    # GOOD UNTIL HERE
    readout_map = pull_curves(canonical_mappings, canonical_replicas, spots_dict)
    # print readout_map
    collector, fails, conservation_object = iterate_through_conditions(readout_map)
    write_out_curves(collector, output_location)
    write_out_curves(fails, total_log)
    # pprint(conservation_object)
    dump(conservation_object, open('cons_obj.dmp', 'w'))

    conservation_object = load(open('cons_obj.dmp', 'r'))

    reduce_table(conservation_object)

    # regress_euploid(conservation_object)
