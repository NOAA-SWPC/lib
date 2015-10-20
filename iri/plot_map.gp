#!/home/palken/usr/bin/gnuplot

set term pngcairo enh col
set out "map.png"

load 'xlonon.cfg'
load 'jet.pal'
load 'lines2.cfg'

set palette maxcolors 0

set yrange [-45:45]
set ytics -45,15,45

set pm3d map interp 20,20
unset key
set label 1 "x10^{12}" right at screen 0.91,0.82
set format cb "%.1f"
set size 0.95,0.95

set cblabel "per m^3"
set ylabel "latitude (degrees)"

set title "O+ density, 12 UT"
splot 'map.dat' us 1:2:($9/10**(12)), 'dat.2000' us 1:2:(0.0) w li lw 4
