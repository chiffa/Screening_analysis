"""
We use the date from the Pavelka, Rancatin and Zhu 2010 Nature paper in order to retrieve the strains that are both
better growing than the wild-type strains and have a lag in their  growth that might signify that there is an adaptation
present.

We use the "significantly-better-growers_2010-09-16.xls" table from the "phenotypic profiling" folder supplied by Jin in
order to determine significantly better growers.

We then scan the data contained in all the subfolders of the phenotypic profiling folder and match the subfolder names
to abbreviations likely to be used to refer to specific drugs

We then write the individual lanes to a file and compute the growth advantage/disadvantage

We then use the data retrieved to peform:
    1) re-combination of replicates into a singe curve+errorbars (TODO)
    2) exponential growth curve regression for each sample (TODO)
    3) Compute the most likely lag + division speed + normalized form (TODO)

"""
__author__ = 'ank@andreikucharavy.com'

import os
import xlrd
import numpy as np
from collections import defaultdict
from csv import reader, writer
from pprint import pprint
from time import sleep

# source_folder = 'L:\\ank\\real-runs\\'
# source_folder = 'L:\\Users\\andrei\\real-runs'
source_folder = 'C:\\Users\\Andrei\\Documents\\real-runs'



# past_mappings = os.path.join(source_folder, 'conditions-to-include-in-clustering_2010-09-16.csv') # used by Pavelka in his clustering
past_mappings = os.path.join(source_folder, 'conditions-to-include-in-clustering_reduced.csv') # Nature paperclustering
# past_mappings = os.path.join(source_folder, 'conditions-to-include-in-clustering_2010-09-16.csv')
better_growers = os.path.join(source_folder, 'significantly-better-growers_2010-09-16.xls')
relative_growth = os.path.join(source_folder, 'meanNormalizedFinalAverageSpotIntensity.csv')
spot_locations = os.path.join(source_folder, 'spotLocation.csv')
output_location = os.path.join(source_folder, 'Re_analysis_ank.csv')
output_folder = os.path.join(source_folder, 'Re_analysis_ank')
name_to_match = 'averageSpotIntensity'

# Conditions name mapping to the
# TODO: remove name collision such as 16C and 16 C
# TODO: inject the initial list with the mappings from Pavelka's mappings
# convert to use rather the names from "conditions to include in clustering"


def build_name2abbreviations(mappings_file):
    storage_dict = defaultdict(list)
    with open(mappings_file, 'rb') as source:
        csv_reader = reader(source)
        csv_reader.next()
        for line in csv_reader:
            if line:
                storage_dict[line[0].split('(')[0]].append(line[2])

    storage_dict = dict([(key, tuple(value)) for key, value in storage_dict.iteritems()])
    return storage_dict


def revert_abbreviation_dict(abbreviation_dict):
    rev_abbreviation_dict = {v: k for k, v in name2abbreviation.items()}
    delkeys = []

    for condition, payload in rev_abbreviation_dict.iteritems():
        if type(condition) == tuple:
            delkeys.append(condition)

    for condition in delkeys:
        for item in condition:
            rev_abbreviation_dict[item] = rev_abbreviation_dict[condition]
        del rev_abbreviation_dict[condition]

    return rev_abbreviation_dict


# recover the better growers from the tables
def get_sig_better_growers(better_growers_table):
    better_pairs = defaultdict(list)

    for sheet in better_growers_table.sheets():
        print 'Sheet:', sheet.name  # just to check we are parsing what we need to parse
        condition = ''
        for row in range(1, sheet.nrows):
            values = []
            if not sheet.cell(row, 0).value:
                values.append(condition)
            else:
                condition = unicode(sheet.cell(row, 0).value)
            for col in range(sheet.ncols):
                v = sheet.cell(row, col).value
                if v is not '':
                    values.append(unicode(v))
            better_pairs[values[0]].append(values[1])

    return better_pairs


def get_better_growers(relative_growth):
    better_pairs = defaultdict(list)
    names = []
    aneuploids = []
    response_table = []
    with open(relative_growth, 'rb') as source:
        csv_reader = reader(source)
        names = csv_reader.next()[1:]
        names = [_name.split('(')[0] for _name in names]
        for line in csv_reader:
            aneuploids.append(line[0])
            response_table.append(line[1:])
    for name in names:
        if name.split('(')[0] not in name2abbreviation.keys():
            print '%s as %s not mapped' % (name, name.split('(')[0])


    aneuploids = np.array(aneuploids)
    response_table = np.array(response_table)
    persistent_shape = response_table.shape
    response_table = np.genfromtxt(response_table.flatten())
    response_table = response_table.reshape(persistent_shape)
    for i, condition in enumerate(names):
        _fltr = response_table[:, i] > 1.01  # TODO: why the heck is this failing?
        better_pairs[condition] += aneuploids[_fltr].tolist()

    better_pairs = dict([(key, list(set(value)))
                         for key, value in better_pairs.iteritems()])
    better_pairs = dict([(key,
                          [_val for _val in value if 'U' not in _val])
                         for key, value in better_pairs.iteritems()])

    return better_pairs


# read location-depended spot locations:
def build_spot_map():
    spot_dict = defaultdict(dict)
    with open(spot_locations, 'rb') as source:
        csv_reader = reader(source)
        for line in csv_reader:
            if line[3]:
                spot_dict[line[3]].update({line[1]: line[2]})
    return spot_dict


def read_line(position, file_path, params=[''], buffer={}):
    if params[0] == (position, file_path) and position in buffer.keys():
        return buffer[position]
    else:
        params[0] = (position, file_path)
        with open(file_path, 'rb') as source:
            csv_reader = reader(source)
            for line in csv_reader:
                buffer[line[0]] = line
                if line[0] == position:
                    return line


# TODO: add replicates grouping
# TODO: add replicates intake from the canonical mappings
# walks the source folder and reads the original data.
def pull_curves(better_pairs, spots_dict):
    locations = []

    stress_specific_locations = defaultdict(list)

    for subfolder in os.walk(source_folder):
        for condition_abbreviation in abbreviation2name.keys():
            parsed_subfolder = subfolder[0].split('\\')
            if len(parsed_subfolder) > 5 \
                    and condition_abbreviation in parsed_subfolder[4] \
                    and parsed_subfolder[-1][-2:] in ['_A', '_B', '_C'] \
                    and os.path.isdir(subfolder[0]):
                file_to_read = ''
                for fle in subfolder[2]:
                    if name_to_match in fle:
                        file_to_read = fle
                condition_name = abbreviation2name[condition_abbreviation]
                plate_type = subfolder[0][-1]
                strains = better_pairs[condition_name]
                # TODO: use plate_type for replicates analysis
                take_field = [spots_dict[strain][plate_type] for strain in strains]
                fpath = subfolder[0]+'\\'+file_to_read
                current_line = read_line('', fpath)[1:]
                current_line.insert(0, '')
                current_line.insert(0, '')
                # current_line.insert(0, '')
                # current_line.insert(0, '')
                sub_locations = []
                sub_locations.append([''])
                # sub_locations.append([fpath,''])
                sub_locations.append(current_line)
                normalized_name = subfolder[0].replace('_A', '')
                normalized_name = normalized_name.replace('_B', '')
                normalized_name = normalized_name.replace('_C', '')
                for i, position in enumerate(take_field):
                    current_line = read_line(position, fpath)[1:]
                    # current_line.insert(0, position)
                    current_line.insert(0, strains[i])
                    # current_line.insert(0, condition_name)
                    current_line.insert(0, parsed_subfolder[4])
                    sub_locations.append(current_line)
                stress_specific_locations[condition_name] += sub_locations
                locations += sub_locations

    return locations, stress_specific_locations

# finally, write out the resulting curves to a destination file
def write_out_curves(locations, out_path):
    with open(out_path, 'wb') as source:
        csv_writer = writer(source)
        for line in locations:
            csv_writer.writerow(line)


if __name__ == "__main__":
    name2abbreviation = {
            '16 C': ('16'),
            'glycerol': ('Glyc'),
            'low glucose': ('Suc','Raf'),
            'high pH': ('hPh'),
            'low pH': ('NaOH', 'lpH'),
            'high salt': ('HS'),
            'rapamycin': ('Rapa'),
            'thiolutin': ('Thio'),
            'fluconazole': ('Fluc'),
            '20 C': ('20C'),
            'bleomycin': ('Bleo'),
            '4-nitroquinoline-N-oxide': ('4NQO'),
            'benomyl': ('Beno')
            }

    name2abbreviation = build_name2abbreviations(past_mappings)
    abbreviation2name = revert_abbreviation_dict(name2abbreviation)

    # good until here

    better_growers_table = xlrd.open_workbook(better_growers)
    # better_pairs = get_sig_better_growers(better_growers_table)
    better_pairs = get_better_growers(relative_growth)
    reference_points = ['controlHaploid', 'controlDiploid', 'controlTriploid']
    better_pairs = dict([(condition, payload + reference_points) for condition, payload in better_pairs.iteritems()])
    spots_dict = build_spot_map()

    locations, stress_specific_locations = pull_curves(better_pairs, spots_dict)
    write_out_curves(locations, output_location)

    for condition, payload in stress_specific_locations.iteritems():
        fname = os.path.join(output_folder, condition + '.csv')
        write_out_curves(payload, fname)

