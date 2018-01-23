#!/usr/bin/env gnuplot

nrow = 3
ncol = 2

load 'multi_default.cfg'

plotheight = 1.5
plotwidth = 5.5
fontsize = ",10"
r = -0.3
l = 0.5
hbuffer = 0.7

load 'multi_defs.cfg'
load 'multi_png.cfg'

outdir = 'output'
file_time = 'window_time.txt'
file_freq = 'window_freq.txt'

load 'grid.cfg'
load 'lines2.cfg'

set point 1

stats file_time

do for [idx=0:int(STATS_blocks - 1)] {

outstr = sprintf('%s/window_%02d.png', outdir, idx)
set out outstr

str = sprintf('Generating plot %s...', outstr)
print str

set multiplot layout nrow,ncol

load 'multi_reset.cfg'

load 'xtimeon.cfg'
set key
set ylabel "{/Symbol \155}A/m^2"

# xtics every 24 hours
set xtics 3600*24

set title "J_r"
plot file_time us 1:2 index idx w lp ti "Original data", \
     file_time us 1:5 index idx w lp ti "Windowed data"

load 'incrow.cfg'

unset key

set title "J_{/Symbol \161}"
plot file_time us 1:3 index idx w lp ti "Original data", \
     file_time us 1:6 index idx w lp ti "Windowed data"

load 'incrow.cfg'

set title "J_{/Symbol \152}"
plot file_time us 1:4 index idx w lp ti "Original data", \
     file_time us 1:7 index idx w lp ti "Windowed data"

load 'inccolumn.cfg'

load 'xtimeoff.cfg'
set xtics auto
set ylabel "Power ({/Symbol \155}A/m^2)"

set title "Amplitude spectrum of J_r"
plot file_freq us 2:3 index idx ls 3 w lp

load 'incrow.cfg'

set title "Amplitude spectrum of J_{/Symbol \161}"
plot file_freq us 2:4 index idx ls 3 w lp

load 'incrow.cfg'

set xlabel "Period (days)"

set title "Amplitude spectrum of J_{/Symbol \152}"
plot file_freq us 2:5 index idx ls 3 w lp

unset xlabel

unset multiplot

}
