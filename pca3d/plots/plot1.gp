#!/usr/bin/env gnuplot

nrow = 2
ncol = 1

load 'multi_default.cfg'

plotwidth = 8.0
fontsize = ",12"
vbuffer = 0.8

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "plot.png"

load 'lines2.cfg'

set multiplot layout nrow,ncol

load 'xtimeon.cfg'
set ylabel "Height-integrated current (A/m)"
set title "Time series of J_{phi} for 110 km altitude, 8 degrees latitude, 150 degrees longitude"
plot 'data_time.txt' us 1:4 w li lt 2 ti "J_{phi}"
load 'xtimeoff.cfg'

load 'incrow.cfg'

set xlabel "Period (hours)"
set ylabel "Power/Frequency (dB/Hz)"
set title "Power Spectral Density of J_{phi} signal"
plot [0:30] 'psd.txt' us 1:2 w li lt 3 ti "PSD of J_{phi}"

unset multiplot
