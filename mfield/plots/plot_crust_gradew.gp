#!/usr/bin/env gnuplot
#
# Plot initial dataset (E/W gradients) used for crustal field modeling

nrow = 2
ncol = 2

load 'multi_default.cfg'

b = 0.5
hbuffer = 1.0

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "crust_gradew.png"

unset key
nskip = 1

load 'xlaton.cfg'
unset xlabel

set multiplot layout nrow,ncol

set yrange [-8:8]

set title "Swarm A and C E/W gradients after removing main and external fields"
set ylabel "dX (nT)"
plot '../output/data1_DX_EW.dat' us 5:(($8-$9) - ($10-$11)) every nskip w p
unset title

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"
set ylabel "dY (nT)"
plot '../output/data1_DY_EW.dat' us 5:(($8-$9) - ($10-$11)) every nskip w p
unset xlabel

load 'inccolumn.cfg'

set ylabel "dZ (nT)"
plot '../output/data1_DZ_EW.dat' us 5:(($8-$9) - ($10-$11)) every nskip w p

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"
set ylabel "dF (nT)"
plot '../output/data1_DF_EW.dat' us 5:(($8-$9) - ($10-$11)) every nskip w p

unset multiplot
