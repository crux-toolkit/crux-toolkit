#!PYTHON
# FILE: compare_sequest_score_spectrum.py
# AUTHOR: CHRIS PARK
# CREATE DATE: NOVEMBER 30 2006

"""
This script compares the results of sequest and CRUX scoring result
"""
import os
import sys
import commands
import random
import re
from pylab import *
from optparse import OptionParser
from parse_sqt_file import Spectrum
from parse_sqt_file import Peptide
from parse_sequest_out_file import SequestObject

#-------------------

def plot_compare_data(result_array, score_type="none", data_type="none"):
    """ plots the results of two compared as a scattered plot """
    prefix = "Score comparison for " + score_type
    file_name = "fig-2-" + data_type + "-" + score_type

    hold(False)
    
    axis('scaled')
    scatter(result_array[0], result_array[1])
    hold(True)
    
    if score_type == "sp":
        #title(prefix + " r = " + r, size=20)
        title(prefix, size=20)
        xlabel("Sequest", size=15)
        ylabel("CRUX", size=15)


        t = range(0, 201)
        y = t
        plot(t, y, color='r')
        
    
        axis([0, 200, 0, 200])
        xticks(size=10)
        yticks(size=10)

    elif score_type == "xcorr":
        title(prefix, size=20)
        xlabel("Sequest", size=15)
        ylabel("CRUX", size=15)
        t = range(0, 4)
        y = t
        plot(t, y, color='r')
        
        xlim((0, 3))
        ylim((0, 3))
    
    hold(False)

    savefig(file_name + ".eps")
    savefig(file_name + ".png")


#-------------------

##needs some work
def plot_histogram_data(data_array, score_type="none"):
    """ plots the results as a histogram """
    prefix = "Score difference for " + score_type
    file_name = "fig-2-difference-" + score_type
    
    #axis('scaled')

    #FIXME Still need a way how to plot all data neatly
    #plot data as histogram
    hold(False)
    hist(data_array, normed=False, bins=10)
    hold(True)
    
    if score_type == "sp":
        title(prefix, size=20)
        # xlim((0, 5))
        xlabel(score_type, size=15)
        ylabel("score difference", size=15)
                
    elif score_type == "xcorr":
        title(prefix, size=20)
        # xlim((0, 5))
        xlabel(score_type, size=15)
        ylabel("score difference", size=15)
        
    savefig(file_name + ".eps")
    savefig(file_name + ".png")
    hold(False)


#-------------------

# Process command line options
usage = "Usage: compare_sequest_score_spectrum <ms2_file> assumes all *.dta files are in current directory "

option_parser = OptionParser(usage)
(options, args) = option_parser.parse_args()

print "number of args: %d\n" % len(args)

if not len(args) == 1:
  print usage
  sys.exit(1)

ms2_file = args[0]
working_dir = "."

# first get all *.dta files to be run and compared
files = os.listdir(working_dir)
pattern = r'\.dta'
dta_file_list = []

#only get list of files with *.dta
for file in files:
    if re.search(pattern, file):
        dta_file_list.append(file)
        
#print files
#print dta_file_list

#score type to compare
score_type_list = ["sp", "xcorr"]

#the result array for each sequest and CRUX run

#for top scored
top_result_array_sp = ([],[])
top_result_array_xcorr = ([],[])

#for randomly selected
random_result_array_sp = ([],[])
random_result_array_xcorr = ([],[])

#Sequest/Crux score difference
diff_sp = []
diff_xcorr = []

totalCount = 0

#seed the random generator with current time
random.seed()

#now iterate all dta files and run sequest
# randomly select one peptide and score against Crux
for dta_file in dta_file_list:
        
    totalCount += 1

    #first get charge from dta_file    
    #get charge state to run comparison
    dta_components = dta_file.rstrip('\n').split('.')
    charge = dta_components[-2]
    scan_number = dta_components[-3]

    
    print "running comparison with charge: %s scan_number: %s" % (charge, scan_number)
    
    #run sequest
    (exit_code, result) = \
                commands.getstatusoutput("./sequest27 " + \
                                         dta_file)
    
    if exit_code == "1":
        print "failed to run sequest"
        sys.exit(1)
        
    #debug
    #print result
    
    #create sequest object
    sequest_object = SequestObject(result)
    
    if sequest_object == None:
        print "failed to create sequest object from sequest output"
        sys.exit(1)
    elif len(sequest_object.peptides) == 0:
        #skip spectrum results with no peptides
        continue

    #peptides to score
    peptide_list = []

    
    # first, randomly select a peptide to score against in crux
    peptide_list.append(random.choice(sequest_object.peptides))
    # second, add first ranked xcorr peptide
    peptide_list.append(sequest_object.peptides[0])

    #now iterate over peptide (random and top peptide)
    for idx, peptide in enumerate(peptide_list):

        #Random peptide?
        if(idx == 0):
            result_array_sp = random_result_array_sp
            result_array_xcorr = random_result_array_xcorr
        else: #top peptide?
            result_array_sp = top_result_array_sp
            result_array_xcorr = top_result_array_xcorr
        
        # for each score type compare
        for score_type in score_type_list:
            #debug
            #print peptide.components

            (exit_code, result_crux) = \
                        commands.getstatusoutput("score_peptide_spectrum " + \
                                                 "--charge " + `charge` + " --score-type " + score_type + " " + \
                                                 peptide.components["sequence"] + " " + \
                                                 `scan_number` + " " + \
                                                 ms2_file 
                                                 ) #add in parameters

            if exit_code == "1":
                print "failed to run score_peptide_spectrum"
                sys.exit(1)

            # print result_crux
        
            result_crux = result_crux.split(': ')

            if score_type == "sp":
                #store both Sequest, Crux scores for the peptide
                result_array_sp[0].append(peptide.components[score_type])
                result_array_sp[1].append(float(result_crux[1]))
                #store the difference in the score values
                diff_sp.append(  abs(float(result_crux[1]) - peptide.components[score_type])  )
            elif score_type == "xcorr":
                result_array_xcorr[0].append(peptide.components[score_type])
                result_array_xcorr[1].append(float(result_crux[1]))
                diff_xcorr.append( abs(float(result_crux[1]) - peptide.components[score_type]) )
        
    if totalCount % 100 == 0:
        print "totalCount: %d" % totalCount
           
#debug
#print result_array_sp
#print result_array_xcorr

# plot the data
print "generating figure"

plot_histogram_data(diff_sp, "sp")
plot_histogram_data(diff_xcorr, "xcorr")
# first plot top ranked peptide comparison
plot_compare_data(top_result_array_sp, "sp", "top")
plot_compare_data(top_result_array_xcorr, "xcorr", "top")

# second, plot randomly picked peptide comparison
plot_compare_data(random_result_array_sp, "sp", "random")
plot_compare_data(random_result_array_xcorr, "xcorr", "random")

# third, plot the score difference of SP, and Xcorr for each peptide between Sequest and Crux
# currently, the data is for both top and random peptides, thus twice much points compared to the other plots
