#!/usr/bin/env python
#
# This is my method of separating data from formatting in a table.
# This script is not intended to be beautiful -- but it should do the
# job.
#
# Run this script to generate table-incremental-optimizations.tex

OUTFILE="table-incremental-optimizations.tex"

tabletext = """

Description
Worm
Yeast
% Change
1
Crux baseline (4/14/09)
5:47:37.9
36:03.5

2
Rewrite including deduplication of peptides; heapify; compressed peptides file; and eliminating seeks.
7:39.0
2:41.6
45.4-fold reduction 
3
Sparse representation of theoretical peaks. [Note that input file parsing is introduced here, and removed later.]

1:50.5
-31.6%
4
Fixed-capacity array for theoretical peaks; better memory management.

1:13.4
-33.3%
5
Linearizing the background subtraction for \\XCorr computation.

0:38.9
-47.2%
6
Active peptide queue and sorted spectra (rolling window join).
0:38.8
0:13.8
-64.5%
7
Eliminate input file parsing. [See text; and compare line 3 above.]
0:25.6

-22.2%
8
Omit calculation with theoretical ions outside spectrometer's range.
0:23.7

-7.4%
9
Sparse difference vector representation.
0:24.7
0:09.0
4.2%
10
Array striping to eliminate one lookup during dot product calculation.
0:20.4

-17.4%
11
Store the vector diffs to disk instead of calculating at runtime.
0:15.5
0:05.9
-24.0%
12
Fixed-point arithmetic.
0:14.7
0:06.0
-5.2%
13
FIFO memory allocator and run-time compiled dot-product code.
0:8.6
0:4.4
-34%
"""

columns = 5
entries = tabletext.splitlines()[1:]
numentries = len(entries)
lines = numentries/columns
assert(lines * columns == numentries)

def entry(line, col):
  return entries[line*columns + col]

def escaped(x):
  return x.replace('%', '\\%')

def MakeLine(fmt, linenum):
  return ' & '.join([fmt % escaped(entry(linenum, col)) for col in range(columns)])

headerline = MakeLine("\\multicolumn{1}{|c|}{%s}", 0)
linesep = """ \\\\
\\hline
"""
otherlines = linesep.join([MakeLine("%s", line) for line in range(1, lines)])

tex = """
\\begin{tabular}{|r|p{0.5\\textwidth}|r|r|p{0.7in}|}
\\hline
%(headerline)s \\\\
\\hline
%(otherlines)s \\\\
\\hline
\\end{tabular}
""" % globals()

open(OUTFILE, "w").write(tex)
