#!/home/palken/usr/bin/gnuplot

set term pngcairo enh col
set out "eval.png"

file = 'eval.txt'

stats file us 1

str = sprintf('maximum eigenvalue = %e', STATS_max)
print str

unset key
set xrange [-0.1:1.1]
unset xtics
xmin = 0.0
xmax = 1.0

set logscale y
set format y "10^{%L}"

set ylabel "normalized eigenvalues"
set title "Spectral Density Matrix Eigenvalues"

set style arrow 1 nohead
plot file us (xmin):(column(1)/STATS_max):(xmax):(column(1)/STATS_max) with vectors arrowstyle 1
