#!/usr/bin/env gnuplot
#
# Plot current stream function with contour lines

nrow = 3
ncol = 1

load 'multi_default.cfg'

plotwidth = 5.0
plotheight = 2.0
b = 0.3

load 'multi_defs.cfg'
load 'multi_png.cfg'

set output "current.png"

file = 'pc_time.txt'

unset key

set xrange [-180:180]
set xtics -180,45,180
#set format x ""

set yrange [-75:75]
set ytics -75,15,75

set cbrange [-10:10]
set cbtics 2

set ylabel "latitude (degrees)"
set cblabel "kA / nT"

# 1 kA spacing between contours
set cntrparam levels inc -10, 1, 10

# Plot contours
do for [idx=1:3] {
  file = sprintf('maps/map%d.dat', idx - 1)
  tfile = sprintf('test%d.dat', idx)
  cfile = sprintf('cont%d.dat', idx)

  set table tfile
  splot file us 1:2:3
  unset table

  set contour base
  unset surface
  set table cfile
  splot file us 1:2:3
  unset table
  unset contour
  set surface
}

set multiplot layout nrow,ncol

do for [idx=1:3] {
  tfile = sprintf('test%d.dat', idx)
  cfile = sprintf('cont%d.dat', idx)

  if (idx == 3) {
    set xlabel "longitude (degrees)"
    set format x "%g"
  } else {
    set format x ""
  }

  load 'moreland.pal'
  set title 'Current stream function, 110 km, PC'.idx
  plot tfile w image, cfile w li lt -1 lw 1.5

  set format x "%g"

  load 'incrow.cfg'
}

unset multiplot
