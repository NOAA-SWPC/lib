#!/usr/bin/env gnuplot
#
# Plot initial dataset (N/S gradients) used for crustal field modeling

nrow = 2
ncol = 2

load 'multi_default.cfg'

b = 0.5
hbuffer = 1.0

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "crust_gradns.png"

unset key
nskip = 1

load 'xlaton.cfg'
unset xlabel

set multiplot layout nrow,ncol

set yrange [-8:8]

set title "Swarm A N/S gradients after removing main and external fields"
set ylabel "dX (nT)"
plot '../output/data0_DX_NS.dat' us 5:(($8-$9) - ($10-$11)) every nskip w p
unset title

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"
set ylabel "dY (nT)"
plot '../output/data0_DY_NS.dat' us 5:(($8-$9) - ($10-$11)) every nskip w p
unset xlabel

load 'inccolumn.cfg'

set ylabel "dZ (nT)"
plot '../output/data0_DZ_NS.dat' us 5:(($8-$9) - ($10-$11)) every nskip w p

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"
set ylabel "dF (nT)"
plot '../output/data0_DF_NS.dat' us 5:(($8-$9) - ($10-$11)) every nskip w p

unset multiplot
