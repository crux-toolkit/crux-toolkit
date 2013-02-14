set output '/dev/null'
set terminal pdf fsize 12
set xzeroaxis
set yzeroaxis
set xtics rotate
set ylabel 'Delta distinct peptide IDs'
set nokey
plot 'plots/weib_triage2.delta.dat' title "Weibull Triage" with points lt 4 pt 4
set output
replot
