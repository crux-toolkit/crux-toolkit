set output "/dev/null"
set terminal png
set title "tide.xcorr"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "tide-xcorr/tide-search.q.txt" using 1:0 title "Tide XCorr" with lines lw 1
replot "tide-xcorr/tide-search.percolator.q.txt" using 1:0 title "Tide XCorr Percolator" with lines lw 1
replot "tide-xcorr/tide-search.q-ranker.q.txt" using 1:0 title "Tide XCorr q-ranker" with lines lw 1
replot "tide-xcorr/tide-search.barista.q.txt" using 1:0 title "Tide XCorr barista" with lines lw 1
set output
replot
