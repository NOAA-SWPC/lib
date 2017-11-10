#!/usr/bin/env gnuplot
#
# Plot current stream function with contour lines

set terminal pngcairo enh col size 1280,960 font ",22"
set output "current.png"

unset key

set xrange [-180:180]
set xtics -180,45,180

set yrange [-75:75]
set ytics -75,15,75

#set cbrange [-100:100]
set cbrange [-60:60]
set cbtics 10

set xlabel "longitude (degrees)"
set ylabel "QD latitude (degrees)"
set cblabel "kA"

set table 'test.dat'
splot 'chi.dat' us 1:2:3
unset table

set contour base
# 10 kA spacing between contours
#set cntrparam levels inc -80, 10, 80
set cntrparam levels inc -80, 5, 80
unset surface
set table 'cont.dat'
splot 'chi.dat' us 1:2:3
unset table

load 'moreland.pal'
set title "Sq current stream function, 110 km"
p 'test.dat' w image, 'cont.dat' w li lt -1 lw 1.5
