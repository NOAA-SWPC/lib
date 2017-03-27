#!/usr/bin/env gnuplot
#
# Plot initial dataset (E/W gradients) used for crustal field modeling

nrow = 3
ncol = 1

load 'multi_default.cfg'

b = 0.5

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
plot '../output/data2_DX_EW.dat' us 5:($8-$9) every nskip w dot
unset title

load 'incrow.cfg'

set ylabel "dY (nT)"
plot '../output/data2_DY_EW.dat' us 5:($8-$9) every nskip w dot

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"
set ylabel "dZ (nT)"
plot '../output/data2_DZ_EW.dat' us 5:($8-$9) every nskip w dot

unset multiplot
