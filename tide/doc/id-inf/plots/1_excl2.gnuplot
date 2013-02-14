set output '/dev/null'
set terminal pdf fsize 12
set xlabel 'Data points available'
set ylabel 'Distinct peptide IDs'
set xtics rotate
set key left top
plot 'plots/e10.dat' notitle with points pt 3 ps 2
replot 'plots/excl2.dat' notitle with points
replot 1199.0*x/7243 notitle lt 1
set output
replot
