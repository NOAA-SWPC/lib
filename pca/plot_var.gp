#!/home/palken/usr/bin/gnuplot

set term pngcairo enh col size 1000,600
set out "variance.png"

file = 'variance_time.txt'

unset key
set xlabel "eigenvalue number"
set ylabel "variance explained"

set point 1
plot [0:70] file w lp lw 4
