#!/usr/bin/env gnuplot
#
# Plot magnetic field response for first few PCs

nrow = 3
ncol = 1

load 'multi_default.cfg'

plotwidth = 5.0
plotheight = 2.0
b = 0.3
hbuffer = 1.5

load 'multi_defs.cfg'
load 'multi_png.cfg'

unset key

set xrange [-180:180]
set xtics -180,45,180

set yrange [-75:75]
set ytics -75,15,75

set ylabel "latitude (degrees)"
set cblabel "dimensionless"

set pm3d map
load 'moreland.pal'

do for [pcidx=1:40] {

  file = sprintf('maps/map%02d.dat', pcidx - 1)
  outfile = sprintf('maps/B%02d.png', pcidx - 1)
  set output outfile

  str = sprintf('Generating plot %s...', outfile)
  print str

  load 'multi_reset.cfg'

  set multiplot layout nrow,ncol

  set format x ""
  unset xlabel

  tstr = sprintf('PC%d', pcidx)

  set title tstr." B_x"
  set cbrange [-10:10]
  splot file us 1:2:4

  load 'incrow.cfg'

  set title tstr." B_y"
  set cbrange [-4:4]
  splot file us 1:2:5

  load 'incrow.cfg'

  set title tstr." B_z"
  set cbrange [-10:10]
  splot file us 1:2:6

  unset multiplot

}
