set output '/dev/null'
set terminal pdf fsize 12
set xzeroaxis
set yzeroaxis
set ylabel 'Delta distinct peptide IDs'
set xlabel 'Delta data points available'
set xtics rotate
set nokey
plot 'plots/chargestates4.delta.dat' title "Charge state exclusion" with points lt 5 pt 5
set output
replot
