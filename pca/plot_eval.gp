#!/usr/bin/env gnuplot

nrow = 1
ncol = 2

load 'multi_default.cfg'

plotheight = 3.0
plotwidth = 5.0
hbuffer = 1.0
r = 0.0
b = 0.2

load 'multi_defs.cfg'

load 'multi_png.cfg'
set out "eval_time.png"

file = 'sval_time.txt'
stats file us 1
sigmax = STATS_max

set multiplot layout nrow,ncol

str = sprintf('maximum singular value = %e', sigmax)
print str

# center of horizontal lines
x = 1.0

xwidth = 1.3
#set xrange [x - xwidth/2 - 0.1:x + xwidth/2 + 0.1]
set xrange [0:2]
set yrange [1e-16:1e1]

unset key
unset xtics

set logscale y
set format y "10^{%L}"

set ylabel "normalized eigenvalues"
set title "covariance matrix spectrum"

set style arrow 1 nohead lw 0.3
plot file us (x - xwidth/2):((column(1)/sigmax)**2):(xwidth):(0.0) with vectors arrowstyle 1

load 'inccolumn.cfg'

unset logscale y
set format y "%g"
set yrange [0:1]
set xtics

set point 1
set xlabel "eigenvalue number"
set ylabel "variance explained"
set title "cumulative variance curve"
plot [0:70] 'variance_time.txt' w lp lw 4 lc rgb "#000000"

unset multiplot
