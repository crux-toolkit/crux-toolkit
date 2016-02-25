set output "/dev/null"
set terminal png
set title "assign-confidence"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "comet/comet.q.txt" using 1:0 title "Comet E-value" with lines lw 1
replot "tide-p-value/tide-search.q.txt" using 1:0 title "Tide p-value" with lines lw 1
replot "tide-xcorr/tide-search.q.txt" using 1:0 title "Tide XCorr" with lines lw 1
set output
replot
