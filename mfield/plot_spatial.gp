#!/home/palken/usr/bin/gnuplot
#
# Plot time/latitude scatter plot of available data

nrow = 3
ncol = 1

load 'multi_default.cfg'

plotwidth = 8.0
b = 0.2

load 'multi_defs.cfg'
load 'multi_png.cfg'
set out "spatial.png"

nskip = 3

set style fill transparent solid 0.2 noborder
set style circle radius 0.004

load 'lines2.cfg'
load 'ylaton.cfg'
#set xrange [2013.5:2016.5]
set xrange [2013.5:*]

satindices = "A B C"

set multiplot layout nrow,ncol

set format x ""
set key tc variable

do for [idx=0:2] {
  file = 'datamap.dat.' . idx
  outstr = sprintf('processing %s...', file)
  print outstr

  if (idx == 2) {
    set format x "%g"
    set xlabel "time"
  }

  sat = word(satindices, idx + 1)
  set title 'Swarm '.sat

  plot file us 1:($6==1 ? $3 : 1/0) every nskip w circles ti "F", \
       file us 1:($9==1 ? $3 : 1/0) every nskip w circles ti "Z", \
       file us 1:($7==1 ? $3 : 1/0) every nskip w circles ti "X", \
       file us 1:($8==1 ? $3 : 1/0) every nskip w circles ti "Y"

  load 'incrow.cfg'

  unset key
}

unset multiplot
