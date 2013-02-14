set output '/dev/null'
set terminal png size 960,720 giant
#set terminal pdf fsize 16 size 9.5,7
#set xrange [-1000:1000]
#set yrange [-300:300]
set xlabel 'Data points available'
set ylabel 'Distinct peptide IDs'
set xtics rotate
set key left top
plot 'excl.dat' title "Parameterized exclusion list" with points lt 1 pt 1
replot 'excl2.dat' title "Basic exclusion list" with points lt 2 pt 2
replot 'e10.dat' title "Actual run (upper right)" with points lt 1 pt 1
replot 1199.0*x/7243 title "Actual IDs per acquired peak" lt 1
set output
replot
