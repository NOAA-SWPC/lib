#!/home/palken/usr/bin/gnuplot

nrow = 1
ncol = 3

load 'multi_default.cfg'

plotheight = 3.0
plotwidth = 5.0
hbuffer = 1.5
r = 0.0
b = 0.5

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "resX.png"

file = 'res.dat'

set pm3d map interp 20,20
unset key

load 'ylaton.cfg'
set xrange [0:360]
set xlabel "longitude (degrees)"
set cbrange [-200:200]
set cblabel "nT"

set multiplot layout nrow,ncol

set title "TIEGCM B_x"
splot file us 1:2:3

load 'inccolumn.cfg'

set title "SH reconstructed B_x"
splot file us 1:2:6

load 'inccolumn.cfg'

set cbrange [*:*]
set title "Residual B_x"
splot file us 1:2:($3-$6)

unset multiplot
