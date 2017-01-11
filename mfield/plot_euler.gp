#!/usr/bin/env gnuplot

nrow = 3
ncol = 1

load 'multi_default.cfg'

plotheight = 1.5
plotwidth = 10.0
vbuffer = 0.8

load 'multi_defs.cfg'

load 'multi_png.cfg'
set out "euler.png"

load 'lines2.cfg'
load 'grid.cfg'

set ylabel "arcseconds"

mylw = 4

set multiplot layout 3,1

set title "Swarm A"
plot 'euler.0' us 2:6 w li lw mylw ti "{/Symbol \141}", \
     'euler.0' us 2:7 w li lw mylw ti "{/Symbol \142}", \
     'euler.0' us 2:8 w li lw mylw ti "{/Symbol \147}"

unset key
load 'incrow.cfg'

set title "Swarm B"
plot 'euler.1' us 2:6 w li lw mylw ti "Alpha", \
     'euler.1' us 2:7 w li lw mylw ti "Beta", \
     'euler.1' us 2:8 w li lw mylw ti "Gamma"

load 'incrow.cfg'

set xlabel "time (years)"
set title "Swarm C"
plot 'euler.2' us 2:6 w li lw mylw ti "Alpha", \
     'euler.2' us 2:7 w li lw mylw ti "Beta", \
     'euler.2' us 2:8 w li lw mylw ti "Gamma"

unset multiplot
