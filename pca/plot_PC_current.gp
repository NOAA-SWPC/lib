#!/usr/bin/env gnuplot
#
# Plot current stream function with contour lines

UT = 12

set term pngcairo enh col size 800,600

unset key
load 'xlonon.cfg'

set yrange [-75:75]
set ytics -75,15,75

set cbrange [-10:10]
set cbtics 2

set ylabel "latitude (degrees)"
set cblabel "kA / nT"

# 1 kA spacing between contours
set cntrparam levels inc -10, 1, 10

tfile = 'test.dat'
cfile = 'cont.dat'

do for [pcidx=1:40] {
  file = sprintf('maps/map_%02dUT_%02d.dat', UT, pcidx)

  outfile = sprintf('maps/current_%02dUT_%02d.png', UT, pcidx)
  set output outfile

  str = sprintf('Generating %s...', outfile)
  print str

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

  load 'moreland.pal'
  set title 'Current stream function, 110 km, '.UT.' UT, PC'.pcidx
  plot tfile w image, cfile w li lt -1 lw 1.5
}
