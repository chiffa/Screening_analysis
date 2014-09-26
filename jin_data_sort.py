__author__ = 'ank'

import os
import xlrd
from collections import defaultdict
from csv import reader, writer
from pprint import pprint

better_growers = 'L:\\ank\\real-runs\\significantly-better-growers_2010-09-16.xls'
spot_locations = 'L:\\ank\\real-runs\\spotLocation.csv'
source_folder = 'L:\\ank\\real-runs'
ret_loc = 'L:\\ank\\real-runs\\Re_analysis_ank.csv'
nm = 'averageSpotIntensity'

wb = xlrd.open_workbook(better_growers)

fw_map ={
            '16 C' : ('16'),
            'glycerol' : ('Glyc'),
            'low glucose' : ('Suc','Raf'),
            'high pH' : ('hPh'),
            'low pH' : ('NaOH', 'lpH'),
            'high salt' : ('HS'),
            'rapamycin' : ('Rapa'),
            'thiolutin' : ('Thio'),
            'fluconazole' : ('Fluc'),
            '20 C' : ('20C'),
            'bleomycin' : ('Bleo'),
            '4-nitroquinoline-N-oxide' : ('4NQO'),
            'benomyl' : ('Beno')
        }

rv_map = {v: k for k, v in fw_map.items()}
delkeys = []
for key, value in rv_map.iteritems():
    if type(key) == tuple:
        delkeys.append(key)

for key in delkeys:
    rv_map[key[0]] = rv_map[key]
    rv_map[key[1]] = rv_map[key]
    del rv_map[key]

better_pairs = defaultdict(list)
for s in wb.sheets():
    print 'Sheet:',s.name
    condition = ''
    for row in range(1, s.nrows):
        values = []
        if not s.cell(row, 0).value:
            values.append(condition)
        else:
            condition = unicode(s.cell(row, 0).value)
        for col in range(s.ncols):
            v = s.cell(row, col).value
            if v is not '':
                values.append(unicode(v))
        better_pairs[values[0]].append(values[1])

spot_dict = defaultdict(dict)
with open(spot_locations, 'rb') as source:
    csv_reader = reader(source)
    for line in csv_reader:
        if line[3]:
            spot_dict[line[3]].update({line[1]: line[2]})

def read_line(position, file_path):
    with open(file_path, 'rb') as source:
        csv_reader = reader(source)
        for line in csv_reader:
            if line[0] == position:
                return line

locations = []

for elt in os.walk(source_folder):
    for cond in rv_map.keys():
        if cond in elt[0] and elt[0][-2:] in ['_A', '_B', '_C'] and os.path.isdir(elt[0]):
            f_elt = ''
            for fle in elt[2]:
                if nm in fle:
                    f_elt = fle
            condition = rv_map[cond]
            plate = elt[0][-1]
            strains = better_pairs[condition]
            take_field = [spot_dict[strain][plate] for strain in strains]
            fpath = elt[0]+'\\'+f_elt
            rl = read_line('', fpath)[1:]
            rl.insert(0, '')
            rl.insert(0, '')
            rl.insert(0, '')
            locations.append([''])
            locations.append([fpath,''])
            locations.append(rl)
            for i, position in enumerate(take_field):
                rl = read_line(position, fpath)[1:]
                rl.insert(0, position)
                rl.insert(0, strains[i])
                rl.insert(0, condition)
                locations.append(rl)


with open(ret_loc, 'wb') as source:
    csv_writer = writer(source)
    for line in locations:
        csv_writer.writerow(line)