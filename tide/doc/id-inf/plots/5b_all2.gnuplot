set output '/dev/null'
set terminal png size 480,720 giant
#set terminal pdf fsize 16 size 4.5,7
#set xrange [-1000:1000]
#set yrange [-300:300]
set xzeroaxis
set yzeroaxis
set xlabel 'Delta data points available'
set ylabel 'Delta distinct peptide IDs'
set xtics rotate
set nokey
plot 'weib_triage2.delta.dat' title "Weibull Triage" with points lt 4 pt 4
replot 'variable2.delta.dat' title "Variable exclusion" with points lt 3 pt 3
replot 'chargestates4.delta.dat' title "Charge state exclusion" with points lt 5 pt 5
set output
replot
