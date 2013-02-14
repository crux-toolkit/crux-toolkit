set output '/dev/null'
set terminal pdf fsize 12
set xrange [4500:7500]
set yrange [750:1200]
set xlabel 'Data points available'
set ylabel 'Distinct peptide IDs'
set xtics rotate
set key left top
plot 'plots/excl2.dat' title "Basic exclusion list" with points lt 2 pt 2
replot 'plots/chargestates4.dat' title "Charge state exclusion" with points lt 5 pt 5
replot 'plots/chargestates3.dat' title "Charge state exclusion, No ID" with points lt 6 pt 6
replot 'plots/e10.dat' notitle with points lt 1 pt 1 ps 2
replot 1199.0*x/7243 notitle lt 1
set output
replot
