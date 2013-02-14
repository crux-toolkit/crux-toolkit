set output '/dev/null'
set terminal png size 480,720 giant
#set terminal pdf fsize 16 size 4.5,7
#set xrange [-1000:1000]
#set yrange [-300:300]
set xlabel 'Data points available'
set ylabel 'Distinct peptide IDs'
set xtics rotate
set key left top
plot 'excl2.dat' title "Basic exclusion list" with points lt 2 pt 2
replot 'weib_triage2.dat' title "Weibull Triage" with points lt 4 pt 4
replot 'variable2.dat' title "Variable exclusion" with points lt 3 pt 3
replot 'chargestates4.dat' title "Charge state exclusion" with points lt 5 pt 5
replot 'e10.dat' title "Actual run (upper right)" with points lt 1 pt 1
replot 1199.0*x/7243 title "Actual IDs per acquired peak" lt 1
set output
replot
