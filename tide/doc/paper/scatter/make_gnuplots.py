#!/usr/bin/env python

gluplot_text = """set output "/dev/null"
set terminal pdf fsize 12
set xlabel "{{Title1}} XCorr"
set ylabel "{{Title2}} XCorr"
set xrange [0:8]
set yrange [0:8]
plot x notitle with lines
replot "scatter/xcorr_xcorr_{{organism}}_{{filename1}}_{{filename2}}.txt" notitle with points
set output
replot
"""

filename_text = "xcorr_xcorr_{{organism}}_{{filename1}}_{{filename2}}.gnuplot"

seqorig = ("seqorig", "SEQUEST 1993")
seqrecent = ("seqrecent", "SEQUEST 2009")
tide = ("tide", "Tide")

pairs = [(seqorig, seqrecent),
         (seqorig, tide),
         (seqrecent, tide)]



import re

def Eval(s, locals):
  return re.sub('{{(.*?)}}', lambda mo:str(eval(mo.group(1), globals(), locals)), s)

for organism in ['yeast', 'worm']:
  for (filename1, Title1), (filename2, Title2) in pairs:
    text = Eval(gluplot_text, globals())
    out_filename = Eval(filename_text, globals())
    open(out_filename, "w").write(text)
