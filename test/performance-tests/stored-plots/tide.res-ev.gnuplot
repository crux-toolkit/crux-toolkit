set output "/dev/null"
set terminal png
set title "tide.res-ev"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "tide-res-ev/tide-search.q.txt" using 1:0 title "Tide res-ev" with lines lw 1
replot "tide-res-ev/tide-search.percolator.q.txt" using 1:0 title "Tide res-ev Percolator" with lines lw 1
set output
replot
