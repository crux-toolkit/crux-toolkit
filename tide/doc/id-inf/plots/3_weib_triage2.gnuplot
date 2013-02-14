set output '/dev/null'
set terminal pdf fsize 12
set xrange [4500:7500]
set yrange [750:1200]
set ylabel 'Distinct peptide IDs'
set xtics rotate
set key left top
plot 'plots/excl2.dat' title "Basic exclusion list" with points lt 2 pt 2
replot 'plots/weib_triage2.dat' title "Weibull Triage" with points lt 4 pt 4
replot 'plots/e10.dat' notitle with points lt 1 pt 1 ps 2
replot 1199.0*x/7243 notitle lt 1
set output
replot
