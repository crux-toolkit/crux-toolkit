set output "/dev/null"
set terminal png
set title "percolator"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "comet/comet.percolator.q.txt" using 1:0 title "Comet Percolator" with lines lw 1
replot "tide-p-value/tide-search.percolator.q.txt" using 1:0 title "Tide p-value Percolator" with lines lw 1
replot "tide-xcorr/tide-search.percolator.q.txt" using 1:0 title "Tide XCorr Percolator" with lines lw 1
set output
replot
