#!/usr/bin/env gnuplot
#
# Plot lat/lon maps of the frequency domain modes contained
# in U_xx.txt data files

nrow = 3
ncol = 2

load 'multi_default.cfg'

vbuffer = 0.5
hbuffer = 1.7
r = -0.5
b = 0.5

load 'multi_defs.cfg'
load 'multi_png.cfg'

outdir = 'output'

set pm3d map interp 0,0
set palette maxcol 0
load 'jet.pal'

load 'ylaton.cfg'

set xrange [0:360]
set format x ""

do for [idx=1:24] {

outstr = sprintf('%s/mode_%02d.png', outdir, idx)
set out outstr

str = sprintf('Generating plot %s...', outstr)
print str

file = sprintf('U_%02d.txt', idx)

set multiplot layout nrow,ncol

load 'multi_reset.cfg'

tstr = sprintf('Real part of radial mode %d', idx)
set title tstr
splot file us 1:2:3

load 'incrow.cfg'

tstr = sprintf('Real part of theta mode %d', idx)
set title tstr
splot file us 1:2:5

load 'incrow.cfg'

set xlabel "longitude (degrees)"
set format x "%g"

tstr = sprintf('Real part of phi mode %d', idx)
set title tstr
splot file us 1:2:7

unset xlabel
set format x ""

load 'inccolumn.cfg'

tstr = sprintf('Imag part of radial mode %d', idx)
set title tstr
splot file us 1:2:4

load 'incrow.cfg'

tstr = sprintf('Imag part of theta mode %d', idx)
set title tstr
splot file us 1:2:6

load 'incrow.cfg'

set xlabel "longitude (degrees)"
set format x "%g"

tstr = sprintf('Imag part of phi mode %d', idx)
set title tstr
splot file us 1:2:8

unset xlabel
set format x ""

unset multiplot

}
