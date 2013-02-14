#!/usr/bin/env python

example_excl="""\
replot 'example.excl.dat' title "Exclusion list" with points lt 2 pt 2
"""

example_variable="""\
replot 'example.variable.dat' title "Variable exclusion" with points lt 3 pt 3
"""

excl="""\
replot 'excl.dat' title "Parameterized exclusion list" with points lt 1 pt 1
"""

excl2="""\
replot 'excl2.dat' title "Basic exclusion list" with points lt 2 pt 2
"""

variable2="""\
replot 'variable2.{"delta." if delta else ""}dat' title "Variable exclusion" with points lt 3 pt 3
"""

weib_triage2="""\
replot 'weib_triage2.{"delta." if delta else ""}dat' title "Weibull Triage" with points lt 4 pt 4
"""

#chargestates2="""\
#replot 'chargestates2.{"delta." if delta else ""}dat' title "Charge state exclusion QR" with points lt 5 pt 5
#"""

chargestates3="""\
replot 'chargestates3.{"delta." if delta else ""}dat' title "Charge state exclusion, No ID" with points lt 1 pt 1
"""

chargestates4="""\
replot 'chargestates4.{"delta." if delta else ""}dat' title "Charge state exclusion" with points lt 5 pt 5
"""

actual_e10="""\
replot 'e10.dat' title "Actual run (upper right)" with points lt 1 pt 1
replot 1199.0*x/7243 title "Actual IDs per acquired peak" lt 1
"""

text="""\
set output '/dev/null'
set terminal png size {960 if name[0]=="1" else 480},720 giant
#set terminal pdf fsize 16 size {9.5 if name[0]=="1" else 4.5},7
#set xrange [-1000:1000]
#set yrange [-300:300]
{"set xzeroaxis\\nset yzeroaxis\\n" if delta else ""}\
set xlabel '{"Delta d" if delta else "D"}ata points available'
set ylabel '{"Delta d" if delta else "D"}istinct peptide IDs'
set xtics rotate
set {"nokey" if delta else "key left top"}
{plots[2:]}\
set output
replot
"""


import re, subprocess

def Eval(s, locals):
  return re.sub('{(.*?)}', lambda mo:str(eval(mo.group(1), globals(), locals)), s)

def Make(name, plots, delta=False):
  t = Eval(Eval(text, locals()), locals())
  gnuplot_prog = "%s.gnuplot" % name
  open(gnuplot_prog, "w").write(t)
  return subprocess.call("gnuplot %s >%s.png" % (gnuplot_prog, name), shell=True)

Make("1_excl_both", excl + excl2 + actual_e10)

Make("1_excl2", excl2 + actual_e10)
Make("2_variable2", excl2 + variable2 + actual_e10)
Make("2b_variable2", variable2, delta=True)
Make("3_weib_triage2", excl2 + weib_triage2 + actual_e10)
Make("3b_weib_triage2", weib_triage2, delta=True)
Make("4_chargestates2", excl2 + chargestates4 + chargestates3 + actual_e10)
Make("4b_chargestates2", chargestates4, delta=True)
Make("5_all2", excl2 + weib_triage2 + variable2 + chargestates4 + actual_e10)
Make("5b_all2", weib_triage2 + variable2 + chargestates4, delta=True)
