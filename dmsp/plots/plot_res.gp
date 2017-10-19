#!/usr/bin/env gnuplot
#
# Plot time series of scalar and Bz residuals before and after calibration

nrow=2
ncol=1

load 'multi_default.cfg'

plotwidth = 10.0
fontsize = ",12"

load 'multi_defs.cfg'
load 'multi_png.cfg'

load 'xtimeon.cfg'
set ylabel "nT"

# DMSP sat number
satnum = "18"
file1 = 'F'.satnum.'_stage1.txt'
file2 = 'F'.satnum.'_stage2.txt'

# line types / colors
LT1 = 5
LT2 = 6

set out 'res_f'.satnum.'_A.png'

set multiplot layout nrow,ncol

set yrange [-200:200]
set key bottom right horizontal tc variable font "Helvetica Bold,18"
set title 'DMSP F-'.satnum.' scalar residuals with CHAOS'
plot file1 us 1:($9-$10) lt LT1 w dot ti "Uncalibrated F"

load 'incrow.cfg'

set yrange [-300:300]
set key bottom right horizontal tc variable font "Helvetica Bold,18"
set title 'DMSP F-'.satnum.' VFM B_3 residuals with CHAOS (projected onto geodetic normal)'
plot file1 us 1:($11-$13) lt LT1 w dot ti "Uncalibrated B_3"

unset multiplot

load 'multi_reset.cfg'
set out 'res_f'.satnum.'_B.png'

set multiplot layout nrow,ncol

idx1=9
idx2=10
set yrange [-200:200]
set title 'DMSP F-'.satnum.' scalar residuals with CHAOS'
plot file1 us 1:($9-$10) lt LT1 w dot ti "Uncalibrated F", \
     file2 us 1:($9-$10) lt LT2 w dot ti "Calibrated F"

load 'incrow.cfg'

idx1=11
idx2=13
set yrange [-300:300]
set title 'DMSP F-'.satnum.' VFM B_3 residuals with CHAOS (projected onto geodetic normal)'
plot file1 us 1:($11-$13) lt LT1 w dot ti "Uncalibrated B_3", \
     file2 us 1:($11-$13) lt LT2 w dot ti "Calibrated B_3"

unset multiplot
