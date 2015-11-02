#!/usr/bin/gnuplot

set term pdfcairo enh col

set out "svplot.pdf"

set pm3d map
unset key

set xrange [-180:180]
set yrange [-90:90]

set xtics (-180,-90,0,90,180)
set ytics (-90,-45,0,45,90)

set xlabel "longitude (degrees)"
set ylabel "latitude (degrees)"
set cblabel "nT"

set cbrange [-150:100]

set palette defined (-150 "dark-blue", -100 "purple", -50 "dark-green", 0 "yellow", 50 "orange", 100 "red")

set multiplot layout 1,2

set colorbox horizontal user origin 0.1,0.15 size 0.8,0.05

set origin 0,0.15
set size 0.5,0.95
set title "WMM 2010 B_z"
splot 'wmm.dat' us 1:2:4

unset colorbox

set origin 0.5,0.15
set size 0.5,0.95
set title "DMSP F16 2008-2010 B_z"
splot 'f16.dat' us 1:2:4

unset multiplot
