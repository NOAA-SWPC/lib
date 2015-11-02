#!/usr/bin/gnuplot

set term pngcairo enh col size 15in,12in font ",10"
set out "dBdt.png"

set style line 1 lt 3

set xrange [2010:*]

set xlabel "time (years)"
set ylabel "nT/year"

set xtics 1

set multiplot layout 3,3 columnsfirst

set key inside top left
set title "KOU dX/dt"
plot 'dBdt_kou.dat' us 1:2 w li ti "DMSP-MAG-1", \
     'kou.txt' us 1:2 w li ls 1 ti "Ground observatory"
unset key

set title "KOU dY/dt"
plot 'dBdt_kou.dat' us 1:3 w li, \
     'kou.txt' us 1:3 w li ls 1

set title "KOU dZ/dt"
plot 'dBdt_kou.dat' us 1:4 w li, \
     'kou.txt' us 1:4 w li ls 1

set title "TAM dX/dt"
plot 'dBdt_tam.dat' us 1:2 w li, \
     'tam.txt' us 1:2 w li ls 1

set title "TAM dY/dt"
plot 'dBdt_tam.dat' us 1:3 w li, \
     'tam.txt' us 1:3 w li ls 1

set title "TAM dZ/dt"
plot 'dBdt_tam.dat' us 1:4 w li, \
     'tam.txt' us 1:4 w li ls 1

set title "MBO dX/dt"
plot 'dBdt_mbo.dat' us 1:2 w li, \
     'mbo.txt' us 1:2 w li ls 1

set title "MBO dY/dt"
plot 'dBdt_mbo.dat' us 1:3 w li, \
     'mbo.txt' us 1:3 w li ls 1

set title "MBO dZ/dt"
plot 'dBdt_mbo.dat' us 1:4 w li, \
     'mbo.txt' us 1:4 w li ls 1

unset multiplot
