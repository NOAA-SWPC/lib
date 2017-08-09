#!/usr/bin/env gnuplot
#
# Plot B_3 residuals vs QD latitude prior to calibration, for both
# CHAOS Bz and CHAOS projected onto geodetic normal

nrow=6
ncol=2

load 'multi_default.cfg'

plotwidth = 6.0
fontsize = ",12"

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "res2.png"

set yrange [-300:300]
load 'xlaton.cfg'
unset xlabel
set key tc variable

set multiplot layout nrow,ncol

idx1=11
idx2=12
set title "VFM B_3 residuals with CHAOS B_z"
do for [i=8:13] {
  istr = sprintf('%02d', i)
  file = 'F15_20'.istr.'_stage1.txt'

  if (i == 13) {
    set xlabel "QD latitude (degrees)"
  }

  plot file us 7:($8 == 1 ? column(idx1)-column(idx2) : 1/0) w dot ti '20'.istr.' Ascending', \
       file us 7:($8 == -1 ? column(idx1)-column(idx2) : 1/0) w dot ti '20'.istr.' Descending'
  unset title
  unset xlabel
  load 'incrow.cfg'
}

load 'inccolumn.cfg'

idx1=11
idx2=13
set title "VFM B_3 residuals with CHAOS projected onto geodetic normal"
do for [i=8:13] {
  istr = sprintf('%02d', i)
  file = 'F15_20'.istr.'_stage1.txt'

  if (i == 13) {
    set xlabel "QD latitude (degrees)"
  }

  plot file us 7:($8 == 1 ? column(idx1)-column(idx2) : 1/0) w dot ti '20'.istr.' Ascending', \
       file us 7:($8 == -1 ? column(idx1)-column(idx2) : 1/0) w dot ti '20'.istr.' Descending'
  unset title
  unset xlabel
  load 'incrow.cfg'
}

unset multiplot
