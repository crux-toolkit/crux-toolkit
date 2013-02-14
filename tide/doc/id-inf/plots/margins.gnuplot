set output '/dev/null'
set terminal pdf fsize 12
set xlabel 'Retention time delta (min)'
set ylabel 'm/z Delta'
set xzeroaxis
set yzeroaxis
set nokey
plot 'plots/margins_out' with dots
set output
replot
