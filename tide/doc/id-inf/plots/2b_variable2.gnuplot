set output '/dev/null'
set terminal pdf fsize 12
set xzeroaxis
set yzeroaxis
set xtics rotate
set ylabel 'Delta distinct peptide IDs'
set nokey
plot 'plots/variable2.delta.dat' title "Variable exclusion" with points lt 3 pt 3
set output
replot
