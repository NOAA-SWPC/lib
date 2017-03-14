#!/usr/bin/env gnuplot
#
# Plot initial dataset used for crustal field modeling

nrow = 3
ncol = 1

load 'multi_default.cfg'

b = 0.5

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "crust_vector.png"

unset key
nskip = 1

load 'xlaton.cfg'
unset xlabel

set multiplot layout nrow,ncol

set yrange [-40:40]

set title "Swarm A vector data after removing main and external fields"
set ylabel "X (nT)"
plot '../output/data0_X.dat' us 5:($7-$8) every nskip w dot
unset title

load 'incrow.cfg'

set ylabel "Y (nT)"
plot '../output/data0_Y.dat' us 5:($7-$8) every nskip w dot

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"
set ylabel "Z (nT)"
plot '../output/data0_Z.dat' us 5:($7-$8) every nskip w dot

unset multiplot
