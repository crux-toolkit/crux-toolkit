set output "/dev/null"
set terminal png
set xlabel "Tide XCorr"
set ylabel "Comet XCorr"
plot x notitle with lines
replot "plots/xcorr.comet.txt" using 4:5 notitle
set output
replot
