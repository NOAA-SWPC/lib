#!/usr/bin/gnuplot
#
# Plot vector residuals from CHAMP showing gravity current spike at -10 QD lat

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
set out "champ_diff_grav.png"

set style line 1 lt 5 lc rgb "red" ps 0.3
set style line 2 lt 5 lc rgb "blue" ps 0.3
set style line 3 lt 5 lc rgb "cyan" ps 0.3

set multiplot layout nrow,ncol columnsfirst

set lmargin at screen left(1,ncol,w,l,r,hbuffer)
set rmargin at screen right(1,ncol,w,l,r,hbuffer)

file1 = 'data/data.dat.2000_2010_15s'
file2 = 'data/data.dat.2000_2010_15s_LT10_14'
file3 = 'data/data.dat.2000_2010_15s_LT15_18'
set yrange [-20:20]
set xrange [-60:60]
set xtics 15
set xtics format ""
set ylabel "along-track difference (nT)"

# First plot

currentplot = 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set key outside top left horizontal tc variable

idx=12
set title "CHAMP B_x 2000-2010, 11 UT, kp <= 2"
plot file1 us 7:($15==1?column(idx):1/0) ls 1 ti "All LT", \
     file2 us 7:($15==1?column(idx):1/0) ls 2 ti "10-14 LT", \
     file3 us 7:($15==1?column(idx):1/0) ls 3 ti "15-18 LT"

unset key

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

idx=13
set title "CHAMP B_y 2000-2010, 11 UT, kp <= 2"
plot file1 us 7:($15==1?column(idx):1/0) ls 1 ti "All LT", \
     file2 us 7:($15==1?column(idx):1/0) ls 2 ti "10-14 LT", \
     file3 us 7:($15==1?column(idx):1/0) ls 3 ti "15-18 LT"

## Last plot

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)
set xlabel "QD latitude (degrees)"
set xtics format "%g"

idx=14
set title "CHAMP B_z 2000-2010, 11 UT, kp <= 2"
plot file1 us 7:($15==1?column(idx):1/0) ls 1 ti "All LT", \
     file2 us 7:($15==1?column(idx):1/0) ls 2 ti "10-14 LT", \
     file3 us 7:($15==1?column(idx):1/0) ls 3 ti "15-18 LT"

unset multiplot
