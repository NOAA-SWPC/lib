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

idx1=9
idx2=10
set key bottom right horizontal tc variable font "Helvetica Bold,18"
set title "Scalar residuals with CHAOS"
plot 'F15_2008_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2008", \
     'F15_2009_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2009", \
     'F15_2010_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2010", \
     'F15_2011_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2011", \
     'F15_2012_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2012", \
     'F15_2013_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2013"

load 'incrow.cfg'

unset key
set yrange [-300:300]

idx1=11
idx2=13
set title "VFM B_3 residuals with CHAOS (projected onto geodetic normal)"
plot 'F15_2008_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2008", \
     'F15_2009_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2009", \
     'F15_2010_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2010", \
     'F15_2011_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2011", \
     'F15_2012_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2012", \
     'F15_2013_stage1.txt' us 1:(column(idx1)-column(idx2)) w dot ti "2013"

load 'incrow.cfg'

set title "VFM H residuals with CHAOS"
plot 'F15_2008_stage1.txt' us 1:(sqrt($9**2 - $11**2)-sqrt($10**2 - $13**2)) w dot ti "2008", \
     'F15_2009_stage1.txt' us 1:(sqrt($9**2 - $11**2)-sqrt($10**2 - $13**2)) w dot ti "2009", \
     'F15_2010_stage1.txt' us 1:(sqrt($9**2 - $11**2)-sqrt($10**2 - $13**2)) w dot ti "2010", \
     'F15_2011_stage1.txt' us 1:(sqrt($9**2 - $11**2)-sqrt($10**2 - $13**2)) w dot ti "2011", \
     'F15_2012_stage1.txt' us 1:(sqrt($9**2 - $11**2)-sqrt($10**2 - $13**2)) w dot ti "2012", \
     'F15_2013_stage1.txt' us 1:(sqrt($9**2 - $11**2)-sqrt($10**2 - $13**2)) w dot ti "2013"

unset multiplot
