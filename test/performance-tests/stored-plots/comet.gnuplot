set output "/dev/null"
set terminal png
set title "comet"
set xlabel "q-value threshold"
set ylabel "Number of accepted PSMs"
set xrange [0:0.1]
set key bottom right
plot "comet/comet.q.txt" using 1:0 title "Comet E-value" with lines lw 1
replot "comet/comet.percolator.q.txt" using 1:0 title "Comet Percolator" with lines lw 1
replot "comet/comet.q-ranker.q.txt" using 1:0 title "Comet q-ranker" with lines lw 1
replot "comet/comet.barista.q.txt" using 1:0 title "Comet barista" with lines lw 1
set output
replot
