#!python
# FILE: sample_spectra_from_ms2.py
# AUTHOR: CHRIS PARK
# CREATE DATE: 8/

"""
samples the spectra on the list of scan number list
print out to stdout
Randomly selects only one Z line thus each spectrum will have only
one charge state
"""
import os.path
import commands
import sys
import random
from optparse import OptionParser

# Process command line options
usage = "Usage: sample_spectra_from_ms2 <ms2_file> <scan number list>"
option_parser = OptionParser(usage)
(options, args) = option_parser.parse_args()

if not len(args) == 2:
  print usage
  sys.exit(1)

ms2_file = args[0]

try:
    scan_num_list = open(args[1], "r")
    
except IOError:
    print "%s\n" % "failed to open file"
    sys.exit(1)

random.seed()

temp_filename = "temp_file"

scan_numbers = []
for line in scan_num_list:
    if line.startswith('S'):
        fields = line.rstrip('\n').split()
        scan_number =  (int)(fields[1].lstrip('0'))
        scan_numbers.append(scan_number)
        

scan_numbers.sort()

for scan_number in scan_numbers:
    
    (exit_code, result) = \
                commands.getstatusoutput("get_ms2_spectrum  " + \
                                         `scan_number` + " " +\
                                         ms2_file + " "+ \
                                         temp_filename
                                         )
    if exit_code == "1":
        print "failed to run"
        sys.exit(1)
        

    try:
        temp_file = open(temp_filename, "r")
        
    except IOError:
        print "%s\n" % "failed to open file"
        sys.exit(1)

    z_num = 0
    z_lines = []
    in_z_line = False
    for t_line  in temp_file:
        if t_line.startswith('File'):
            continue
        t_line = t_line.rstrip('\n')
        if t_line.startswith('Z'):
            in_z_line = True
            z_num += 1
            z_lines.append(t_line)
            continue
        elif in_z_line:
            in_z_line = False
            select_z_line = random.choice(z_lines)
            print select_z_line
            
        print t_line
                
        
    (exit_code, result) = \
                commands.getstatusoutput("rm -rf " + temp_filename)
        
                                             
