#!/usr/bin/env python

# for plotting ideas:
# ~/tide/fdr2.py
# http://matplotlib.sourceforge.net/plot_directive/mpl_examples/pylab_examples/bar_stacked.py

# DATA; Times in seconds for Worm 10k

profile_tidezero = """
                        Tide-v0
Dot Product             221.3
Preprocess Spectrum     115.7
Theoretical Peaks       101.5
Other                   20.7     
"""

profiles_linearize = """
                        Before Linearizing                   After Linearizing
Dot Product             8.074                                8.558
Preprocess Spectrum     35.966                               1.167
Theoretical Peaks       27.158                               27.23
Other                   2.202                                1.945
"""

profiles_fivefold = """
                        Before 5x sparser    After 5x sparser    After Storing Diffs
Dot Product             12.087               7.41                6.0064
Preprocess Spectrum     1.3035               1.976               1.877
Theoretical Peaks       6.399                1.976               1.3139
Other                   3.9105               6.175               9.5727
Get Sparse Peaks                             7.163
"""

profile_final = """
                        Tide
Dot Product             1.40
Preprocess Spectrum     1.58
Theoretical Peaks       0.894
Other                   2.56
"""

import re, matplotlib.pyplot as plot

def GetFields(line, col_start_end):
  return [f.strip() for f in [line[start:end] for start, end in col_start_end]]

def GetTable(t):
  two_blank_splitter = re.compile(r'\S(?:\S|(?:\s(?!\s|$)))+') # split string by minimum of 2-blanks
  t = t.split('\n')[1:-1]
  col_seps = [m.start() for m in two_blank_splitter.finditer(t[0])]
  col_start_end = zip([0] + col_seps, col_seps + [100000])
  fields = [GetFields(line, col_start_end) for line in t]
  col_labels = fields[0][1:]
  row_labels = [f[0] for f in fields[1:]]
  vals = [[float('0%s' % entry) for entry in row[1:]] for row in fields[1:]]
  return (col_labels, row_labels, vals)

def nearest_round(x, r):
  r += 0.0
  return int(x/r + 0.5) * r

def Plot(plt, table, title, ylabel = True):
  (col_labels, row_labels, vals) = table
  num_cols = len(col_labels)
  width = 0.1       # the width of the bars: can also be len(x) sequence
  ind = [i/(num_cols+1.0) for i in range(1, num_cols+1)]    # the x locations
  colors = ['blue', 'red', 'green', 'yellow', 'purple']

  bottom = [0] * num_cols
  p = []
  for vrow, color in zip(vals, colors):
    p.append(plt.bar([i-width/2.0 for i in ind], vrow, width, color = color, bottom = bottom))
    bottom = [i+j for i, j in zip(bottom, vrow)]

  if ylabel: plt.ylabel('Time (s)')
  plt.title(title, horizontalalignment = 'left', x = 0)
  plt.xlim(0,1)
  col_labels = [label.replace(' ', '\n') for label in col_labels]
  plt.xticks(ind, col_labels, multialignment='center', size='small')
  largest = max(bottom)
  plt.ylim(0, largest)
  y_tick_divisions = 5
  y_tick_size = largest / (y_tick_divisions+0.0)
  rounding = 0.1 if largest < 100 else 1
  plt.yticks([nearest_round(y_tick_size * i, rounding) for i in range(y_tick_divisions + 1)], size='small')
  return ([q[0] for q in p], row_labels) # for legend


plot.subplot(221)
Plot(plot, GetTable(profile_tidezero), '(a)')
plot.subplot(222)
Plot(plot, GetTable(profiles_linearize), '(b)', ylabel = False)
plot.subplot(223)
rects, labels = Plot(plot, GetTable(profiles_fivefold), '(c)')
plot.subplot(224)
Plot(plot, GetTable(profile_final), '(d)', ylabel = False)
plot.subplots_adjust(top=0.95, left=0.1, right=0.75, hspace=0.3)
labels = [label.replace(' ', '\n') for label in labels]
labels[-1] = 'Get Sparse\nPeaks'
plot.figlegend(rects, labels, 'center right')
plot.savefig('profiles.png')
