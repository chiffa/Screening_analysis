"""
We use the date from the Pavelka, Rancatin and Zhu 2010 Nature paper in order to retrieve the strains that are both
better growing than the wild-type strains and have a lag in their  growth that might signify that there is an adaptation
present.

We use the "significantly-better-growers_2010-09-16.xls" table from the "phenotypic profiling" folder supplied by Jin in
order to determine significantly better growers.

We then scan the data contained in all the subfolders of the phenotypic profiling folder and match the subfolder names
to abbreviations likely to be used to refer to specific drugs


"""
__author__ = 'ank@stowers.org'

import os
import xlrd
import numpy as np
from collections import defaultdict
from csv import reader, writer
from pprint import pprint

better_growers = 'L:\\ank\\real-runs\\significantly-better-growers_2010-09-16.xls'
relative_growth = 'L:\\ank\\real-runs\\meanNormalizedFinalAverageSpotIntensity.csv'
spot_locations = 'L:\\ank\\real-runs\\spotLocation.csv'
source_folder = 'L:\\ank\\real-runs'
output_location = 'L:\\ank\\real-runs\\Re_analysis_ank.csv'
output_folder = 'L:\\ank\\real-runs\\Re_analysis_ank'
name_to_match = 'averageSpotIntensity'

# Conditions name mapping to the
# TODO: remove name collision such as 16C and 16 C
# TODO: inject the initial list with the mappings from Pavelka's mappings
# convert to use rather the names from "conditions to include in clustering"
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

# perform the reverse mapping of the names
abbreviation2name = {v: k for k, v in name2abbreviation.items()}

# some keys that are tuples and need to be unpacked to be properly available
delkeys = []

for condition, payload in abbreviation2name.iteritems():
    if type(condition) == tuple:
        delkeys.append(condition)

for condition in delkeys:
    abbreviation2name[condition[0]] = abbreviation2name[condition]
    abbreviation2name[condition[1]] = abbreviation2name[condition]
    del abbreviation2name[condition]
#############################################################################
# end of block

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
    table = []
    with open(relative_growth, 'rb') as source:
        csv_reader = reader(source)
        names = csv_reader.next()[1:]
        for line in csv_reader:
            aneuploids.append(line[0])
            table.append(line[1:])
    for name in names:
        if name.split('(')[0] not in abbreviation2name.keys():
            print '%s as %s not mapped' % (name, name.split('(')[0])

    aneuploids = np.array(aneuploids)
    table = np.array(table)
    for i, condition in enumerate(names):
        fltr = table[i,:] > 1
        print condition, aneuploids[fltr]
        better_pairs[condition] += aneuploids[fltr].tolist()

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

# walk the source folder
# TODO: add matching to the strains that were used for clustering in the original paper
# TODO: add the link to all the other better growers
# TODO: add replicates grouping
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
                plate = subfolder[0][-1]
                strains = better_pairs[condition_name]
                take_field = [spots_dict[strain][plate] for strain in strains]
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
                normalized_name = normalized_name.replace('_A', '')
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
    better_growers_table = xlrd.open_workbook(better_growers)
    better_pairs = get_sig_better_growers(better_growers_table)
    get_better_growers(relative_growth)
    reference_points = ['controlHaploid', 'controlDiploid', 'controlTriploid']
    better_pairs = dict([(condition, payload + reference_points) for condition, payload in better_pairs.iteritems()])
    spots_dict = build_spot_map()
    pprint(dict(spots_dict))
    locations, stress_specific_locations = pull_curves(better_pairs, spots_dict)
    write_out_curves(locations, output_location)

    for condition, payload in stress_specific_locations.iteritems():
        fname = os.path.join(output_folder, condition + '.csv')
        write_out_curves(payload, fname)

