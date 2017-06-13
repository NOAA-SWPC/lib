#!/usr/bin/env gnuplot
#
# Plot time series of scalar, Bz, and H residuals prior to calibration

nrow=3
ncol=1

load 'multi_default.cfg'

plotwidth = 10.0
fontsize = ",12"

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "res1.png"

load 'xtimeon.cfg'
set ylabel "nT"
set yrange [-200:200]

hypot(x,y) = sqrt(x*x + y*y)

set multiplot layout nrow,ncol

idx1=3
idx2=7
set key bottom right horizontal tc variable font "Helvetica Bold,18"
set title "Scalar residuals with CHAOS"
plot 'dat08' us 1:(column(idx1)-column(idx2)) w dot ti "2008", \
     'dat09' us 1:(column(idx1)-column(idx2)) w dot ti "2009", \
     'dat10' us 1:(column(idx1)-column(idx2)) w dot ti "2010", \
     'dat11' us 1:(column(idx1)-column(idx2)) w dot ti "2011", \
     'dat12' us 1:(column(idx1)-column(idx2)) w dot ti "2012", \
     'dat13' us 1:(column(idx1)-column(idx2)) w dot ti "2013"

load 'incrow.cfg'

unset key
set yrange [-300:300]

idx1=6
idx2=10
set title "VFM B_3 residuals with CHAOS (projected onto geodetic normal)"
plot 'dat08' us 1:(column(idx1)-column(idx2)) w dot ti "2008", \
     'dat09' us 1:(column(idx1)-column(idx2)) w dot ti "2009", \
     'dat10' us 1:(column(idx1)-column(idx2)) w dot ti "2010", \
     'dat11' us 1:(column(idx1)-column(idx2)) w dot ti "2011", \
     'dat12' us 1:(column(idx1)-column(idx2)) w dot ti "2012", \
     'dat13' us 1:(column(idx1)-column(idx2)) w dot ti "2013"

load 'incrow.cfg'

set title "VFM H residuals with CHAOS"
plot 'dat08' us 1:(hypot($4,$5)-hypot($8,$9)) w dot ti "2008", \
     'dat09' us 1:(hypot($4,$5)-hypot($8,$9)) w dot ti "2009", \
     'dat10' us 1:(hypot($4,$5)-hypot($8,$9)) w dot ti "2010", \
     'dat11' us 1:(hypot($4,$5)-hypot($8,$9)) w dot ti "2011", \
     'dat12' us 1:(hypot($4,$5)-hypot($8,$9)) w dot ti "2012", \
     'dat13' us 1:(hypot($4,$5)-hypot($8,$9)) w dot ti "2013"

unset multiplot
