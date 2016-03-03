set output "/dev/null"
set terminal png
set title "q-ranker"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "comet/comet.q-ranker.q.txt" using 1:0 title "Comet q-ranker" with lines lw 1
replot "tide-xcorr/tide-search.q-ranker.q.txt" using 1:0 title "Tide XCorr q-ranker" with lines lw 1
replot "tide-p-value/tide-search.q-ranker.q.txt" using 1:0 title "Tide p-value q-ranker" with lines lw 1
set output
replot
