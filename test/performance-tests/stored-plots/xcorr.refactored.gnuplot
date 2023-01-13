set output "/dev/null"
set terminal png
set xlabel "Tide XCorr"
set ylabel "Refactored XCorr"
plot x notitle with lines
replot "plots/xcorr.refactored.txt" using 4:5 notitle
set output
replot
