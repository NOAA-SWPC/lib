#!/usr/bin/env gnuplot

nrow = 3
ncol = 2

load 'multi_default.cfg'

plotheight = 1.5
plotwidth = 5.5
fontsize = ",10"

load 'multi_defs.cfg'
load 'multi_png.cfg'

outdir = 'output'
file_time = 'window_time.txt'
file_freq = 'window_freq.txt'

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
set xtics 3600*8
set key
set ylabel "{/Symbol \155}A/m^2"

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

plot file_freq us 

unset multiplot

}
