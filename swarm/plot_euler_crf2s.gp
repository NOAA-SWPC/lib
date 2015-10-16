#!/usr/bin/gnuplot

### n: change this parameter to equal the number of data sets to be plotted
nrow = 3
ncol = 1

# Height of each plot in pixels
plotheight = 150.0

# Width of each plot in pixels
plotwidth = 500.0

# Vertical buffer space between each plot in pixels
vbuffer = 20

# Horizontal buffer space between each plot in pixels
hbuffer = 120

# t: top margin in pixels
t = 50.0

# b: key height in pixels (bottom margin)
b = 70.0

# h: height of total plot in pixels
h = (plotheight+vbuffer)*nrow + t + b

# l: left margin in pixels
l = 100.0

# r: right margin in pixels
r = -50.0

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
set out "crf2s.png"

set ylabel "angle (degrees)"
set xrange [2013.9:2015]

set style line 1 lt 1 lw 2
set style line 2 lt 2 lw 2
set style line 3 lc rgb "blue" lw 2
set style line 4 lt 4 lw 2
set style line 5 lt 5 lw 2
set style line 6 lc rgb "orange" lw 2
set style line 7 lc rgb "brown" lw 2
set style line 8 lc rgb "black" lw 2

set xtics 0.2 format "" nomirror

set multiplot layout (nrow+0),ncol columnsfirst

set lmargin at screen left(1,ncol,w,l,r,hbuffer)
set rmargin at screen right(1,ncol,w,l,r,hbuffer)

unset key

# First plot

set yrange [-4000:4000]

currentplot = 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \141} (arcseconds)"
plot 'crf2s.dat.A' us 1:5 w li ls 1 ti "A", \
     'crf2s.dat.B' us 1:5 w li ls 2 ti "B", \
     'crf2s.dat.C' us 1:5 w li ls 3 ti "C"

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \142} (arcseconds)"
plot 'crf2s.dat.A' us 1:6 w li ls 1 ti "A", \
     'crf2s.dat.B' us 1:6 w li ls 2 ti "B", \
     'crf2s.dat.C' us 1:6 w li ls 3 ti "C"

## Last plot

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)
set xlabel "time (years)"
set xtics format "%g"

set key at screen 0.5,0.05 horizontal center
set ylabel "{/Symbol \104}{/Symbol \147} (arcseconds)"
plot 'crf2s.dat.A' us 1:7 w li ls 1 ti "A", \
     'crf2s.dat.B' us 1:7 w li ls 2 ti "B", \
     'crf2s.dat.C' us 1:7 w li ls 3 ti "C"

unset key
unset xlabel

unset multiplot
