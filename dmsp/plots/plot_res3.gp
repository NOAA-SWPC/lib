#!/usr/bin/env gnuplot
#
# Plot scalar residuals vs QD latitude

nrow=6
ncol=1

load 'multi_default.cfg'

plotwidth = 6.0
fontsize = ",12"

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "res3.png"

set yrange [-300:300]
load 'xlaton.cfg'
unset xlabel
set key tc variable

set multiplot layout nrow,ncol

idx1=3
idx2=7
set title "DMSP F-15 scalar residuals against CHAOS"
do for [i=8:13] {
  istr = sprintf('%02d', i)
  file = 'dat'.istr

  if (i == 13) {
    set xlabel "QD latitude (degrees)"
  }

  plot file us 2:($12 == 1 ? column(idx1)-column(idx2) : 1/0) w dot ti '20'.istr.' Ascending', \
       file us 2:($12 == -1 ? column(idx1)-column(idx2) : 1/0) w dot ti '20'.istr.' Descending'
  unset title
  unset xlabel
  load 'incrow.cfg'
}

unset multiplot
