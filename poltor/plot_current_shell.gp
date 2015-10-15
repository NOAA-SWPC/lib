#!/usr/bin/env gnuplot
#
# Plot current stream function on higher shell (350km) with contour lines

set terminal pngcairo enh col
set output "current_shell.png"

unset key

set xrange [-180:180]
set xtics -180,45,180

set yrange [-75:75]
set ytics -75,15,75

set cbrange [-100:100]

set xlabel "longitude (degrees)"
set ylabel "QD latitude (degrees)"
set cblabel "kA"

set table 'test.dat'
splot 'chi.dat' us 1:2:4
unset table

set contour base
# 10 kA spacing between contours
set cntrparam levels inc -80, 10, 80
unset surface
set table 'cont.dat'
splot 'chi.dat' us 1:2:4
unset table

set pm3d map interp 20,20
set palette @MATLAB
#splot 'chi.dat' us 1:2:4
p 'test.dat' w image, 'cont.dat' w li lt -1 lw 1.5
