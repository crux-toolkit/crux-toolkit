#!/usr/bin/env python

# Fields:
# Program ('crux', 'tide', 'sequest', etc.)
# Organism ('yeast', 'worm') or with 2 phosphorylations, ('yeast_phos2', 'worm_phos2')
# Precursor mass tolerance window (0.25, 3.00)
# Digestion ('full', 'partial')
# Average spectra/second over all replicate runs
# Number of spectra
# Average run time over all replicate runs
# Number of replicate runs done (3 or 5)
# List of runtimes for each replicate run

results = [
('crux', 'worm', 3.00, 'full', 17.26, 10000, 579.43, 3, [596.31, 571.07, 570.92]),
('crux', 'worm', 3.00, 'partial', 1.01, 100, 99.03, 5, [114.80, 95.34, 95.09, 94.81, 95.09]),
('crux', 'worm', 0.25, 'full', 67.40, 10000, 148.36, 5, [148.80, 148.39, 148.16, 148.32, 148.14]),
('crux', 'worm', 0.25, 'partial', 9.95, 100, 10.05, 5, [10.11, 10.03, 10.03, 10.05, 10.03]),
('crux', 'yeast', 3.00, 'full', 42.03, 10000, 237.93, 3, [239.98, 237.14, 236.67]),
('crux', 'yeast', 3.00, 'partial', 3.70, 100, 27.05, 5, [27.25, 27.13, 26.95, 26.95, 26.96]),
('crux', 'yeast', 0.25, 'full', 99.90, 10000, 100.10, 5, [100.13, 99.92, 99.91, 99.90, 100.62]),
('crux', 'yeast', 0.25, 'partial', 29.02, 100, 3.45, 5, [3.51, 3.42, 3.44, 3.44, 3.42]),
('tide', 'worm', 3.00, 'full', 1554.00, 10000, 6.43, 5, [6.46, 6.40, 6.51, 6.42, 6.38]),
('tide', 'worm', 0.25, 'full', 1979.57, 10000, 5.05, 5, [5.12, 5.02, 5.04, 5.03, 5.05]),
('tide', 'worm', 3.00, 'partial', 77.23, 10000, 129.48, 5, [128.89, 130.80, 128.91, 128.47, 130.32]),
('tide', 'worm', 0.25, 'partial', 412.93, 10000, 24.22, 5, [24.27, 24.18, 24.16, 24.25, 24.23]),
('tide', 'yeast', 3.00, 'full', 2737.63, 10000, 3.65, 5, [3.68, 3.62, 3.68, 3.66, 3.63]),
('tide', 'yeast', 0.25, 'full', 3243.80, 10000, 3.08, 5, [3.09, 3.08, 3.09, 3.08, 3.07]),
('tide', 'yeast', 3.00, 'partial', 536.29, 10000, 18.65, 5, [18.92, 18.60, 18.60, 18.87, 18.25]),
('tide', 'yeast', 0.25, 'partial', 1198.44, 10000, 8.34, 5, [8.35, 8.27, 8.50, 8.26, 8.35]),
('tide', 'worm_phos2', 3.00, 'full', 318, 10000, 31.438, 3, [27.062, 27.521, 39.732]),
('tide', 'yeast_phos2', 3.00, 'full', 960, 10000, 10.413, 3, [8.662, 8.589, 13.989]),
('sequest_orig', 'worm', 0.25, 'full', 0.77, 100, 130.00, 5, [130.42, 130.54, 129.67, 129.69, 129.68]),
('sequest_orig', 'worm', 3.00, 'full', 0.63, 100, 158.23, 5, [157.65, 157.73, 158.29, 157.81, 159.69]),
('sequest_orig', 'yeast', 0.25, 'full', 2.83, 100, 35.34, 5, [35.38, 35.35, 35.23, 35.36, 35.36]),
('sequest_orig', 'yeast', 3.00, 'full', 1.89, 100, 52.94, 5, [55.80, 52.31, 52.36, 51.99, 52.23]),
('sequest_orig', 'worm_phos2', 3.00, 'full', 0.503, 100, 198.99, 3, [198.986, 199.0, 198.983]),
('sequest_orig', 'yeast_phos2', 3.00, 'full', 1.21, 100, 82.633, 3, [82.442, 83.544, 81.913]),
('omssa', 'yeast', 0.25, 'full', 206.41, 10000, 48.45, 5, [48.545, 48.485, 48.432, 48.453, 48.317]),
('omssa', 'worm', 0.25, 'full', 84.64, 10000, 118.15, 5, [118.054, 118.322, 118.011, 118.148, 118.193]),
('omssa', 'worm', 3.0, 'full', 44.70, 10000, 223.73, 5, [223.516, 224.853, 223.556, 223.394, 223.334]),
('omssa', 'yeast', 3.0, 'full', 124.00, 10000, 80.64, 5, [80.402, 81.754, 80.334, 80.514, 80.208]),
('omssa', 'worm', 3.0, 'partial', 0.76, 10000, 13107.47, 4, [13107.3, 13118.3, 13103.6, 13100.7]),
('omssa', 'yeast', 3.0, 'partial', 2.97, 10000, 3363.74, 3, [3366.29, 3362.16, 3362.76]),
('omssa', 'worm', 0.25, 'partial', 7.21, 10000, 1386.20, 3, [1385.04, 1386.32, 1387.23]),
('omssa', 'yeast', 0.25, 'partial', 26.60, 10000, 375.94, 3, [376.052, 376.051, 375.732]),
('xtandem', 'yeast', 0.25, 'full', 399.27, 10000, 25.046, 1, [25.046]),
('xtandem', 'yeast', 0.25, 'partial', 109.28, 10000, 91.512, 1, [91.512]),
('xtandem', 'yeast', 3.00, 'full', 138.17, 10000, 72.374, 1, [72.374]),
('xtandem', 'yeast', 3.00, 'partial', 12.39, 10000, 807.389, 1, [807.389]),
('xtandem', 'worm', 0.25, 'full', 198.83, 10000, 50.293, 1, [50.293]),
('xtandem', 'worm', 0.25, 'partial', 33.33, 10000, 300.031, 1, [300.031]),
('xtandem', 'worm', 3.00, 'full', 46.29, 10000, 216.037, 1, [216.037]),
('xtandem', 'worm', 3.00, 'partial', 3.19, 10000, 3139.47, 1, [3139.47]),
('sequest_indexed', 'worm', 3.00, 'full', 8.77, 10000, 1139.76, 3, [1125.2, 1172.22, 1121.87]),
('sequest_indexed', 'yeast', 3.00, 'full', 24.75, 10000, 404.12, 3, [403.817, 404.114, 404.43]),
('sequest_indexed', 'worm', 0.25, 'full', 55.43, 10000, 180.41, 3, [175.635, 179.038, 186.557]),
('sequest_indexed', 'yeast', 0.25, 'full', 104.4, 10000, 95.78, 3, [95.734, 95.763, 95.849]),
('sequest_indexed', 'worm', 3.00, 'partial', 1.62, 10000, 6171.54, 3, [6166.49, 6170.78, 6177.35]),
('sequest_indexed', 'yeast', 3.00, 'partial', 6.16, 10000, 1623.75, 3, [1631.12, 1622.98, 1617.16]),
('sequest_indexed', 'worm', 0.25, 'partial', 15.18, 10000, 658.74, 3, [648.345, 647.519, 680.352]),
('sequest_indexed', 'yeast', 0.25, 'partial', 42.8, 10000, 233.64, 3, [229.151, 228.452, 243.324]),
('sequest_indexed', 'worm_phos2', 3.00, 'full', 1.22, 100, 82.02, 3, [81.409, 81.568, 83.084]),
('sequest_indexed', 'yeast_phos2', 3.00, 'full', 3.91, 100, 25.584, 3, [25.553, 25.243, 25.956])
# ('sequest_recent', 'worm', 0.25, 'full', 1.31, 100, 76.43, 5, [74.29, 85.00, 74.27, 74.31, 74.27]),
# ('sequest_recent', 'worm', 0.25, 'partial', 2.32, 100, 43.09, 5, [43.09, 43.07, 43.04, 43.17, 43.08]),
# ('sequest_recent', 'worm', 3.00, 'full', 1.29, 100, 77.23, 5, [76.05, 76.03, 81.95, 76.04, 76.06]),
# ('sequest_recent', 'worm', 3.00, 'partial', 0.89, 100, 111.98, 5, [112.16, 113.27, 110.45, 110.55, 113.47]),
# ('sequest_recent', 'yeast', 0.25, 'full', 4.90, 100, 20.41, 5, [20.41, 20.40, 20.41, 20.40, 20.42]),
# ('sequest_recent', 'yeast', 0.25, 'partial', 8.58, 100, 11.66, 5, [11.66, 11.70, 11.65, 11.65, 11.65]),
# ('sequest_recent', 'yeast', 3.00, 'full', 4.75, 100, 21.06, 5, [21.01, 21.10, 20.95, 20.97, 21.28]),
# ('sequest_recent', 'yeast', 3.00, 'partial', 3.34, 100, 29.91, 5, [29.88, 29.93, 29.92, 29.89, 29.90]),
# ('sequest_recent', 'worm_phos2', 3.00, 'full', 0.81, 100, 123.405, 3, [123.321, 123.579, 123.314]),
# ('sequest_recent', 'yeast_phos2', 3.00, 'full', 2.78, 100, 35.963, 3, [36.568, 36.097, 35.224]),
]

f1 = [('sequest_orig',    'SEQUEST'),
      ('crux',            'Crux'),
      ('omssa',           'OMSSA'),
      ('sequest_indexed', 'Indexed SEQUEST'),
      ('xtandem',         'X!Tandem'),
      ('tide',            'Tide')]

f1_phos = [('sequest_orig',    'SEQUEST'),
           ('sequest_indexed', 'Indexed SEQUEST'),
           ('tide',            'Tide')]

f2 = [('partial', 'Semi-tryptic digest'),
      ('full',    'Fully tryptic digest')]

f3 = [('worm',  'Worm'),
      ('yeast', 'Yeast')]

f3_phos = [('worm_phos2',  'Worm'),
           ('yeast_phos2', 'Yeast')]

f4 = [(3,    '\xb13.0 Da'),
      (0.25, '\xb10.25 Da')]

def roun(x):
  if (x >= 100):
    return int(round(x))
  if (x >= 10):
    return round(x*10)/10.0
  return x

d = {}
for z1, z2, z3, z4, z5, z6, z7, z8, z9 in results:
  d[(z1, z2, z3, z4)] = roun(z5)

import itertools

def ColHeaders(*fields):
  headers = []
  for field in range(len(fields)):
    line = ['']
    last = None
    for i in itertools.product(*fields):
      dummy, f = i[field]
      line.append(f if f != last else '')
      last = f
    headers.append(line)
  return headers


def RowData(lookup_order, row_field, *col_fields):
  col_fields = [[x for x, y in f] for f in list(col_fields)]
  lines = []
  for row_f, row_f_header in row_field:
    line = [row_f_header]
    for i in itertools.product(*col_fields):
      j = [row_f] + list(i)
      k = [j[k[0]] if type(k) == type([]) else k for k in lookup_order]
      line.append(d.get(tuple(k),''))
    lines.append(line)
  return lines


import sys, csv

def WriteTable(filename, field_list):
  w = csv.writer(open(filename, "w"))
  for i in field_list:
    w.writerow(i)
  

WriteTable("timing_results_main.csv", ColHeaders(f2, f3, f4) + RowData([[0], [2], [3], [1]], f1, f2, f3, f4))
WriteTable("timing_results_mods.csv", ColHeaders(f3_phos) + RowData([[0], [1], 3.0, 'full'], f1_phos, f3_phos))



# progs = [('sequest_orig', 'SEQUEST', '7/1993'),
#          ('sequest_recent', 'SEQUEST', '11/2009'),
#          ('crux', '', 'Crux'),
#          ('tide', '', 'Tide')
#          ]
# numprogs = len(progs)
# 
# result_dict = {}
# for prog, organism, mz_tol, partial_or_full, specs_per_sec, spectra, time, num_runs, raw in results:
#   result_dict[(prog, organism, mz_tol, partial_or_full)] = (specs_per_sec, spectra)
# 
# from math import log10
# 
# def Entry(prog, organism, mz_tol, partial_or_full):
#   entry = result_dict.get((prog[0], organism, mz_tol, partial_or_full), None)
#   if not entry:
#     return ""
#   specs_per_sec, spectra = entry
#   digits = 2-int(log10(specs_per_sec))
#   pattern = "$%%0.%df$" % (digits if digits >= 0 else 0)
#   return pattern % (specs_per_sec) # "^*" if spectra == 100 else ""
# 
# def GetLines():
#   result = ""
#   for organism in ['worm', 'yeast']:
#     organism_header = "\\hline\n\\multirow{4}{*}{%s}" % organism
#     for digest in ['partial', 'full']:
#       digest_header = "\\multirow{2}{*}{%s}" % digest
#       for MZ_TOLERANCE in [0.25, 3.00]:
#         result += ("%s & %s & %0.2f & %s \\\\\n"
#                    % (organism_header, digest_header, MZ_TOLERANCE, 
#                       " & ".join([Entry(p, organism, MZ_TOLERANCE, digest)  for p in progs])))
#         organism_header = ""
#         digest_header = ""
#   return result
# 
# # double brackets [[expression]] will be replaced with the value of the expression
# 
# text = """
# \\begin{tabular}{ccc|[["c"*numprogs]]}
#           &        &   mass    & [[" & ".join([formalname1 for label, formalname1, formalname2 in progs])]] \\\\
# benchmark & digest & tolerance & [[" & ".join([formalname2 for label, formalname1, formalname2 in progs])]] \\\\
# 
# [[GetLines()]]
# 
# \\end{tabular}
# """
# 
# import re
# print re.sub(r'[[][[](.*?)[]][]]', lambda m: "%s" % eval(m.group(1), globals()), text) 
