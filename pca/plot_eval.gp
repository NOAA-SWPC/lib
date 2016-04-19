#!/home/palken/usr/bin/gnuplot

set term pngcairo enh col
set out "eval.png"

#file = 'eval.txt'
file = 'eval_time.txt'

stats file us 1

str = sprintf('maximum eigenvalue = %e', STATS_max)
print str

# center of horizontal lines
x = 2.0

xwidth = 1.0
set xrange [x - xwidth/2 - 0.1:x + xwidth/2 + 0.1]
set yrange [*:1e5]

unset key
unset xtics

set logscale y
set format y "10^{%L}"

set ylabel "normalized eigenvalues"
set title "Spectral Density Matrix Eigenvalues"

set style arrow 1 nohead lw 0.4
plot file us (x - xwidth/2):(column(1)/STATS_max):(xwidth):(0.0) with vectors arrowstyle 1
