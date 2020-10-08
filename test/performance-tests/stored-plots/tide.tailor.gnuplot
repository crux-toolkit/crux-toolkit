set output "/dev/null"
set terminal png
set title "tide.tailor"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "tide-tailor/tide-search.q.txt" using 1:0 title "Tide Tailor " with lines lw 1
replot "tide-tailor/tide-search.percolator.q.txt" using 1:0 title "Tide Tailor Percolator" with lines lw 1
replot "tide-tailor/tide-search.q-ranker.q.txt" using 1:0 title "Tide Tailor q-ranker" with lines lw 1
replot "tide-tailor/tide-search.barista.q.txt" using 1:0 title "Tide Tailor barista" with lines lw 1
set output
replot
