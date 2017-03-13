#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 16 March 2009
import sys
import os
import math

usage = """USAGE: make-qq-plot.py <p-values> <root>

Compare a given set of p-values to the uniform distribution by
creating a QQ plot with log-log axes.  The program outputs three
files: a gnuplot script (<root>.gnuplot), the data to be plotted
(<root>.txt) and the plot itself (<root>.png).  Note that the stored
values are downsampled to avoid having too many points in the plot.


Options:
  --no-log-scale
  --column-header <string>  Header of column from which to get p-values.
  --minus-natural-log       Input values are negative log base e.
  --format png|eps          (default=png)
  --fontsize <int>          (only effective with "-format eps")
  --title <title>

If the p-value file is specified as "-", then the program reads from
standard input.

"""

###############################################################################
# Find a given word in a tab-delimited string of words.
# Return the index.
def findWord(header, word):

  words = header.split("\t")
  for index in range(0, len(words)):
    if (words[index] == word):
      return(index)
  sys.stderr.write("Can't find %s in %s.\n" % (word, header))
  sys.exit(1)

###############################################################################
# MAIN
###############################################################################

# Set default values.
log_scale = 1
column_header = ""
log_values = 0
file_format = "png"
font_size = 24
title = ""

# Parse the command line.
sys.argv = sys.argv[1:]
while (len(sys.argv) > 2):
  next_arg = sys.argv[0]
  sys.argv = sys.argv[1:]
  if (next_arg == "--no-log-scale"):
    log_scale = 0
  elif (next_arg == "--column-header"):
    column_header = sys.argv[0]
    sys.argv = sys.argv[1:]
  elif (next_arg == "--minus-natural-log"):
    log_values = 1
  elif (next_arg == "--format"):
    file_format = sys.argv[0]
    sys.argv = sys.argv[1:]
  elif (next_arg == "--fontsize"):
    font_size = int(sys.argv[0])
    sys.argv = sys.argv[1:]
  elif (next_arg == "--title"):
    title = sys.argv[0]
    sys.argv = sys.argv[1:]
  else:
    sys.stderr.write("Invalid option (%s).\n" % next_arg)
    sys.exit(1)
if (len(sys.argv) != 2):
  sys.stderr.write(usage)
  sys.exit(1)
pvalue_filename = sys.argv[0]
fileroot = sys.argv[1]

# Open the file for reading.
if (pvalue_filename == "-"):
  pvalue_file = sys.stdin
else:
  pvalue_file = open(pvalue_filename, "r")

# If a header string was specified, find the relevant column.
if (column_header != ""):
  header = pvalue_file.readline().rstrip()
  column_index = findWord(header, column_header)
  sys.stderr.write("Reading p-values from column %d.\n" % column_index)
else:
  column_index = 0

# Read the p-values from the specified column.
pvalues = []
numZeroes = 0
for line in pvalue_file:
  line = line.rstrip()
  words = line.split("\t")

  # Skip comment lines.
  if (line[0] == "#"):
    continue

  # Crash if the line is too short.
  if (len(words) <= column_index):
    sys.stderr.write("Too few columns (%d < %d).\n%s\n" 
                     % (len(words), column_index, line))
    sys.exit(1)

  # Skip NaNs.
  if ((words[column_index] == "NaN") or
      (words[column_index] == "nan")):
    continue

  pvalue = float(words[column_index])
  if (log_values):
    pvalue = math.exp(-1.0 * pvalue)

  # Count zero p-values.
  if (pvalue == 0):
    numZeroes += 1

  # Store this p-value.
  pvalues.append(pvalue)

pvalue_file.close()
num_pvalues = len(pvalues)
if (numZeroes != 0):
  sys.stderr.write("Warning: Found %d zero p-values.\n" % numZeroes)
sys.stderr.write("Read %d p-values from %s.\n" % (num_pvalues, 
                                                  pvalue_filename))

# Sort the values.
pvalues.sort()

# Open the data file.
data_filename = "%s.txt" % fileroot
data_file = open(data_filename, "w")
sys.stderr.write("Creating %s.\n" % data_filename)

# We will only print with this density along the x-axis.
if (log_scale):
  increment = 0.01
else:
  increment = 0.001
current_value = 0

# Print the values to a file.
rank = 1.0
num_printed = 0
for pvalue in pvalues:

  if (log_scale):
    new_value = math.log(rank / num_pvalues)
  else:
    new_value = rank / num_pvalues

  if (current_value == 0) or (new_value >= current_value + increment):
    data_file.write("%g\t%g\n" % (rank / num_pvalues, pvalue))
    current_value = new_value
    num_printed += 1

  rank += 1.0
data_file.close()
sys.stderr.write("Printed %d p-values.\n" % num_printed)

# Find the first non-zero p-value.
for index in range(0, len(pvalues)):
  min_pvalue = pvalues[index]
  if (min_pvalue != 0):
    break

# Set the range.
sys.stderr.write("Minimum p-value=%g\n" % min_pvalue)
if (1.0 / num_pvalues < min_pvalue):
  min_pvalue = 1.0 / num_pvalues
  sys.stderr.write("Minimum rank p-value=%g\n" % min_pvalue)
if (min_pvalue == 0):
  min_value = "1e-10"
else:
  min_value = "1e%d" % (int(math.log(min_pvalue, 10.0)) - 1)
sys.stderr.write("Minimum x-axis value=%s\n" % min_value)

# Open the gnuplot file.
gnuplot_filename = "%s.gnuplot" % fileroot
gnuplot_file = open(gnuplot_filename, "w")
sys.stderr.write("Creating %s.\n" % gnuplot_filename)

# Print the gnuplot file.
gnuplot_file.write("set output '/dev/null'\n")
if (file_format == "png"):
  gnuplot_file.write("set terminal png\n")
elif (file_format == "eps"):
  gnuplot_file.write("set terminal postscript eps %s\n" % font_size)
else:
  sys.stderr.write("Invalid file format (%s).\n" % file_format)
  sys.exit(1)
gnuplot_file.write("set xlabel 'Rank p-value'\n")
gnuplot_file.write("set ylabel 'Calculated p-value'\n")
gnuplot_file.write("set xrange [%s:1]\n" % min_value)
gnuplot_file.write("set yrange [%s:1]\n" % min_value)
if (log_scale):
  gnuplot_file.write("set logscale xy\n")
if (title != ""):
  gnuplot_file.write("set title '%s'\n" % title)
gnuplot_file.write("plot x notitle with lines lt 1\n")
gnuplot_file.write("replot 0.5*x notitle with lines lt 2\n")
gnuplot_file.write("replot 2.0*x notitle with lines lt 2\n")
gnuplot_file.write("replot '%s' notitle with points\n" % data_filename)
gnuplot_file.write("set output\n")
gnuplot_file.write("replot\n")
gnuplot_file.close()

# Make the image.
sys.stderr.write("Creating %s.%s.\n" % (fileroot, file_format))
os.system("gnuplot %s > %s.%s" % (gnuplot_filename, fileroot, file_format))

