set output "/dev/null"
set terminal png
set title "tide.combined"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "tide-combined/tide-combined.q.txt" using 1:0 title "Tide combined p-value" with lines lw 1
replot "tide-combined/tide-combined.percolator.q.txt" using 1:0 title "Tide combined p-value Percolator" with lines lw 1
set output
replot
