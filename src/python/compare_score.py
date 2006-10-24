#!/usr/bin/env python2.4
# FILE: compare_score.py
# AUTHOR: CHRIS PARK
# CREATE DATE: 23 OCTOBER 2006

"""
this script compares the scores from the .sqt file and my file
"""

import os.path
import commands
import sys
import parse_sqt_file
from optparse import OptionParser

# Process command line options
usage = "Usage: compare_score <score_type> <sqt_file> <ms2_file>"
option_parser = OptionParser(usage)
(options, args) = option_parser.parse_args()

if not len(args) == 3:
  print usage
  sys.exit(1)

#set sp_score type
score_type = args[0]

# add more score
if not score_type == "sp":
    print usage
    sys.exit(1)

#check if sqt file and ms2 can be accesseed

sqt_object = SqtObject(args[1])

if sqt_object == None:
    sys.exit(1);


result_array = []

for working_spectrum in sqt_object.spectrums:
    scanNum = working_spectrum.fields["id"]
    for working_peptide in working_spectrum.peptides:
        (exit_code, result) = \
                    commands.getstatusoutput('score_peptide_spectrum ') #add in parameters
        result_array.append( (#real result, result))

plot #
