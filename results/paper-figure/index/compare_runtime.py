#!PYTHON
# FILE: compare_runtime.py
# AUTHOR: CHRIS PARK
# CREATE DATE: 7/30/2007

"""
This script compares the runtime between Crux and Sequest
"""
import os
import sys
import commands
from scipy import *
from pylab import *
from optparse import OptionParser

#-------------------

def plot_compare_data(crux_array, sequest_array, mass_windows, number_of_spectrum, score_type="xcorr"):
    """compares runtime for each scoring method """

    print "generating figures"
    
    xlabel("mass window (Da)", size=15)
    ylabel("Runtime for " + `number_of_spectrum` + " spectra (s)", size=15)
    
    plot(mass_windows, crux_array, label="Crux")
    plot(mass_windows, sequest_array, label="SEQUEST")
    legend()
    
    savefig("indexing" + ".eps")
    savefig("indexing" + ".png")


#_________________________________________
#
# By, Christopher Park
# This script compares the runtime between of Sequest and CRUX
#_________________________________________


# Process command line options
usage = "Usage: compare_runtime <ms2 file> <fasta_file>"
option_parser = OptionParser(usage)
(options, args) = option_parser.parse_args()

print "number of args: %d\n" % len(args)

if not len(args) == 2:
  print usage
  sys.exit(1)

ms2_file = args[0]
fasta_file = args[1]

#mass windows to test runtime
mass_windows = ["0.1", "1", "3"]

#runtime result arrays for each method
sequest_results = []
crux_results = []

number_of_spectrum = 0

charge_list = [1,2,3]

#first run Sequest with varying mass windows
for window in mass_windows:

    # 1, run Sequest
    seq_time = 0.0
    for charge in charge_list:
        (exit_code, result) = \
                    commands.getstatusoutput("time -p ./sequest27 " + \
                                             "-Psequest.params" + "_" + window + \
                                             " " + "*." + "$charge.dta"
                                             )
        #debug
        print result
        #print exit_code
    
        if exit_code == "1":
            print "%s %s" % ("failed to run Sequest on mass window:", window)
            sys.exit(1)
        else:
            #now parse the runtime from the result output
            result = result.split('\n')
            for line in result:
                #get user time
                #FIXME is it user or real?
                if line.startswith('real '):
                    fields = line. rstrip('\n').split()
                    seq_time += float(fields[1])
    sequest_results.append(seq_time)

                
    # 2, now run Crux
    (exit_code, result) = \
                commands.getstatusoutput("time -p search_spectra " + \
                                         "--mass-window " + window + " "  \
                                         "--parameter-file " + "crux_parameter "
                                         + `ms2_file` + " " \
                                         + `fasta_file`
                                         )
    #debug
    print result
    #print exit_code
    
    if exit_code == "1":
        print "%s %s" % ("failed to run Crux on mass window:", window)
        sys.exit(1)
    else:
        #now parse the runtime from the result output
        result = result.split('\n')
        for line in result:
            #get user time
            if line.startswith('user '):
                fields = line.rstrip('\n').split()
                crux_results.append(float(fields[1]))
            elif line.startswith('# SPECTRUM SCAN'):
                number_of_spectrum += 1

#Debug
#print crux_results, sequest_results, mass_windows

#now plot the results
plot_compare_data(crux_results, sequest_results, [ float(i) for i in mass_windows], number_of_spectrum)
