set output '/dev/null'
set terminal pdf fsize 12
set xrange [0:100]
set yrange [0:100]
set xlabel 'e10 retention time (min)'
set ylabel 'e12 retention time (min)'
plot ((0.994061 * x) + -0.452867) title 'y = 0.994x + -0.453'
replot 'plots/e10-vs-e12.dups' notitle with points
set output
replot
