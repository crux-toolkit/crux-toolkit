#!PYTHON
# FILE: compare_score.py
# AUTHOR: CHRIS PARK
# CREATE DATE: 23 OCTOBER 2006

"""
this script compares the scores from the .sqt file and my file
"""
from pylab import *
import os.path
import commands
import sys
from optparse import OptionParser
from parse_sqt_file import SqtObject
from parse_sqt_file import Spectrum
from parse_sqt_file import Peptide

#-------------------

def plot_compare_data(result_array, score_type="sp"):
    """ plots the results of two compared as a scattered plot """
    prefix = "Sequest-vs-CRUX-for-" + score_type
    
    title(prefix, size=20)
    xlabel("Sequest", size=15)
    ylabel("CRUX", size=15)

    # plot each point
    #for (sequest, crux) in result_array:
    scatter(result_array[0], result_array[1])

    # plot y=x
    t = range(0, 300)
    y = t
    plot(t, y, color='r')
    
    #legend(loc='lower right')
    axis([0, 300, 0, 300])
    xticks(size=10)
    yticks(size=10)

    savefig(prefix + ".eps")
    savefig(prefix + ".png")

#-------------------

def plot_compare_rank(result_array, score_type="sp"):
    prefix = "Sequest_Xcore_rank-vs-CRUX-for-" + score_type
    title(prefix, size=20)
    xlabel("Sequest Xcore rank", size=15)
    ylabel("CRUX", size=15)

    # plot each point
    #for (sequest, crux) in result_array:
    scatter(result_array[0], result_array[1])

    # plot y=x
    #t = range(0, 300)
    #y = t
    #plot(t, y, color='r')
    
    #legend(loc='lower right')
    axis([0, 300, 0, 300])
    xticks(size=10)
    yticks(size=10)

    savefig(prefix + ".eps")
    savefig(prefix + ".png")



    
#-------------------

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
sqt_file = args[1]
ms2_file = args[2]

sqt_object = SqtObject(sqt_file)

if sqt_object == None:
    sys.exit(1)


result_array = ([],[])
result_array2 = ([],[])
totalCount = 0

for working_spectrum in sqt_object.spectrums:
    if totalCount >= 3000:
        break
    scanNum = working_spectrum.fields["id"]
    charge = working_spectrum.fields["charge"]
    for working_peptide in working_spectrum.peptides:
        (exit_code, result) = \
                    commands.getstatusoutput("score_peptide_spectrum " + \
                                             "--charge " + `charge` + " --score_type " + score_type + " " + \
                                             working_peptide.components["sequence"] + " " + \
                                             `scanNum` + " " + \
                                             ms2_file 
                                             ) #add in parameters
        if exit_code == "1":
            print "failed to run score_peptide_spectrum"
            sys.exit(1)
        #(real result, result) store in result array
        #for line in result:
        #    if line.startswith('I'):
        #        print line
        #        continue
        #    elif line.startswith('S'):
        result = result.split(': ')
        #if working_peptide.components["xcore_rank"] < 10:
        result_array[0].append(working_peptide.components[score_type])
        result_array[1].append(float(result[1]))
        #    result_array2[0].append(working_peptide.components["xcore_rank"])
        #    result_array2[1].append(float(result[2]))
        totalCount += 1
            #if working_peptide.components[score_type] >= 210:
            #print "Sequest: %.2f, CRUX: %.2f, sequence: %s" % (working_peptide.components[score_type], float(result[1]),working_peptide.components["sequence"])
        if totalCount % 100 == 0:
            print "totalCount: %d" % totalCount
            #break
        
        
#plot the data
plot_compare_data(result_array, score_type)
#plot_compare_rank(result_array2)

#print each pair of score
#(result1, result2) = result_array
#n = 0
#for stuff in result1:
#    print "%.1f %.1f" % (stuff, result2[n])
#    n += 1
