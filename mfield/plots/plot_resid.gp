#!/usr/bin/env gnuplot
#

nrow = 4
ncol = 1

load 'multi_default.cfg'

plotheight = 2.0
plotwidth = 10.0
vbuffer = 0.8

load 'multi_defs.cfg'
load 'multi_png.cfg'

satnum = '0'
iternum = '1'
set out 'res'.satnum.'_iter'.iternum.'.png'

unset key
nskip = 10

load 'xlaton.cfg'

unset xlabel
set ylabel "nT"

set multiplot layout nrow,ncol

file = '../output/res'.satnum.'_X_iter'.iternum.'.dat'
set title "X residuals"
plot file us 5:($8-$9-$10) every nskip w p

load 'incrow.cfg'

file = '../output/res'.satnum.'_Y_iter'.iternum.'.dat'
set title "Y residuals"
plot file us 5:($8-$9-$10) every nskip w p

load 'incrow.cfg'

file = '../output/res'.satnum.'_Z_iter'.iternum.'.dat'
set title "Z residuals"
plot file us 5:($8-$9-$10) every nskip w p

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"

file = '../output/res'.satnum.'_F_iter'.iternum.'.dat'
set title "F residuals"
plot file us 5:($8-$9) every nskip w p

unset multiplot
