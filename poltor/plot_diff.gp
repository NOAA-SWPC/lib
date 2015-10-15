#!/usr/bin/gnuplot

### n: change this parameter to equal the number of data sets to be plotted
nrow = 3
ncol = 1

# Height of each plot in pixels
plotheight = 200.0

# Width of each plot in pixels
plotwidth = 1000.0

# Vertical buffer space between each plot in pixels
vbuffer = 50

# Horizontal buffer space between each plot in pixels
hbuffer = 80

# t: top margin in pixels
t = 75.0

# b: key height in pixels (bottom margin)
b = 50.0

# h: height of total plot in pixels
h = (plotheight+vbuffer)*nrow + t + b

# l: left margin in pixels
l = 100.0

# r: right margin in pixels
r = 0.0

# w: width of total plot in pixels
w = ncol*plotwidth + l + r

### define functions to help set top/bottom margins
# h-t-b               = plotting area vertical height (excluding top/bottom margins)
# (h-t-b)/n           = height of each plot including buffer = plotheight + vbuffer
# (h-t-b)*(i-1)/n     = vertical position of plot i relative to 0
# t + (h-t-b)*(i-1)/n = vertical position of plot i relative to t
top(i,n,h,t,b,buf) = 1.0 - (t+(h-t-b)*(i-1)/n)/h
bot(i,n,h,t,b,buf) = 1.0 - (t+(h-t-b)*i/n-buf)/h

left(j,n,w,l,r,buf) = (l + (w-l-r)/n*(j-1))/w
right(j,n,w,l,r,buf) = (l + (w-l-r)/n*j-buf)/w

set term pngcairo enh col size w,h font ",10"
set out "champ_diff.png"

set style line 1 lt 1 lw 2
set style line 2 lt 2 lw 2
set style line 3 lc rgb "blue" lw 2
set style line 4 lt 4 lw 2
set style line 5 lt 5 lw 2
set style line 6 lc rgb "orange" lw 2
set style line 7 lc rgb "brown" lw 2
set style line 8 lc rgb "black" lw 2

set multiplot layout nrow,ncol columnsfirst

set lmargin at screen left(1,ncol,w,l,r,hbuffer)
set rmargin at screen right(1,ncol,w,l,r,hbuffer)

file = 'data/data.dat.2000_2010_15s'
set yrange [-80:80]
set xrange [-90:90]
set xtics 15
set xtics format ""

qdidx = 8

# First plot

currentplot = 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set key outside top left horizontal tc variable

set title "CHAMP B_x 2000-2010, 11 UT, kp <= 2"
plot file us qdidx:($16==1&$17==0?$10:1/0) w dot ti "B", \
     file us qdidx:($16==1&$17==0?$13:1/0) w dot ti "dB" lt 3

unset key

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "CHAMP B_y 2000-2010, 11 UT, kp <= 2"
plot file us qdidx:($16==1&$17==0?$11:1/0) w dot, \
     file us qdidx:($16==1&$17==0?$14:1/0) w dot lt 3

## Last plot

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)
set xlabel "QD latitude (degrees)"
set xtics format "%g"

set title "CHAMP B_z 2000-2010, 11 UT, kp <= 2"
plot file us qdidx:($16==1&$17==0?$12:1/0) w dot, \
     file us qdidx:($16==1&$17==0?$15:1/0) w dot lt 3

unset multiplot
