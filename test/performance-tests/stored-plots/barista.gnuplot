set output "/dev/null"
set terminal png
set title "barista"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "comet/comet.barista.q.txt" using 1:0 title "Comet barista" with lines lw 1
replot "tide-p-value/tide-search.barista.q.txt" using 1:0 title "Tide p-value barista" with lines lw 1
replot "tide-xcorr/tide-search.barista.q.txt" using 1:0 title "Tide XCorr barista" with lines lw 1
set output
replot
