set output "/dev/null"
set terminal png
set title "tide.p-value"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "tide-p-value/tide-search.q.txt" using 1:0 title "Tide p-value" with lines lw 1
replot "tide-p-value/tide-search.percolator.q.txt" using 1:0 title "Tide p-value Percolator" with lines lw 1
replot "tide-p-value/tide-search.q-ranker.q.txt" using 1:0 title "Tide p-value q-ranker" with lines lw 1
replot "tide-p-value/tide-search.barista.q.txt" using 1:0 title "Tide p-value barista" with lines lw 1
set output
replot
