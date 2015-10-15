#!/usr/bin/gnuplot

### n: change this parameter to equal the number of data sets to be plotted
nrow = 11
ncol = 3

# Height of each plot in pixels
plotheight = 150.0

# Width of each plot in pixels
plotwidth = 400.0

# Vertical buffer space between each plot in pixels
vbuffer = 70

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
set out "gd_champ.png"

set ylabel "residual (nT)"
set xrange [-50:50]
set yrange [-20:20]
set xtics -50, 10, 50

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

unset key

# First plot

currentplot = 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

# QD idx
qdidx = 10

# X index
idx = 17

set title "2000 CHAMP X residuals"
plot 'data/champ_2000.dat' us qdidx:idx w dot, \
     'data/champ_2000_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2001 CHAMP X residuals"
plot 'data/champ_2001.dat' us qdidx:idx w dot, \
     'data/champ_2001_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2002 CHAMP X residuals"
plot 'data/champ_2002.dat' us qdidx:idx w dot, \
     'data/champ_2002_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2003 CHAMP X residuals"
plot 'data/champ_2003.dat' us qdidx:idx w dot, \
     'data/champ_2003_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2004 CHAMP X residuals"
plot 'data/champ_2004.dat' us qdidx:idx w dot, \
     'data/champ_2004_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2005 CHAMP X residuals"
plot 'data/champ_2005.dat' us qdidx:idx w dot, \
     'data/champ_2005_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2006 CHAMP X residuals"
plot 'data/champ_2006.dat' us qdidx:idx w dot, \
     'data/champ_2006_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2007 CHAMP X residuals"
plot 'data/champ_2007.dat' us qdidx:idx w dot, \
     'data/champ_2007_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2008 CHAMP X residuals"
plot 'data/champ_2008.dat' us qdidx:idx w dot, \
     'data/champ_2008_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2009 CHAMP X residuals"
plot 'data/champ_2009.dat' us qdidx:idx w dot, \
     'data/champ_2009_X.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set xlabel "QD latitude (degrees)"
set title "2010 CHAMP X residuals"
plot 'data/champ_2010.dat' us qdidx:idx w dot, \
     'data/champ_2010_X.dat' us 1:2 w li ls 3
unset xlabel

# Y index
idx = 18

set lmargin at screen left(2,ncol,w,l,r,hbuffer)
set rmargin at screen right(2,ncol,w,l,r,hbuffer)

currentplot = 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2000 CHAMP Y residuals"
plot 'data/champ_2000.dat' us qdidx:idx w dot, \
     'data/champ_2000_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2001 CHAMP Y residuals"
plot 'data/champ_2001.dat' us qdidx:idx w dot, \
     'data/champ_2001_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2002 CHAMP Y residuals"
plot 'data/champ_2002.dat' us qdidx:idx w dot, \
     'data/champ_2002_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2003 CHAMP Y residuals"
plot 'data/champ_2003.dat' us qdidx:idx w dot, \
     'data/champ_2003_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2004 CHAMP Y residuals"
plot 'data/champ_2004.dat' us qdidx:idx w dot, \
     'data/champ_2004_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2005 CHAMP Y residuals"
plot 'data/champ_2005.dat' us qdidx:idx w dot, \
     'data/champ_2005_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2006 CHAMP Y residuals"
plot 'data/champ_2006.dat' us qdidx:idx w dot, \
     'data/champ_2006_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2007 CHAMP Y residuals"
plot 'data/champ_2007.dat' us qdidx:idx w dot, \
     'data/champ_2007_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2008 CHAMP Y residuals"
plot 'data/champ_2008.dat' us qdidx:idx w dot, \
     'data/champ_2008_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2009 CHAMP Y residuals"
plot 'data/champ_2009.dat' us qdidx:idx w dot, \
     'data/champ_2009_Y.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set xlabel "QD latitude (degrees)"
set title "2010 CHAMP Y residuals"
plot 'data/champ_2010.dat' us qdidx:idx w dot, \
     'data/champ_2010_Y.dat' us 1:2 w li ls 3
unset xlabel

# Z index
idx = 19

set lmargin at screen left(3,ncol,w,l,r,hbuffer)
set rmargin at screen right(3,ncol,w,l,r,hbuffer)

currentplot = 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2000 CHAMP Z residuals"
plot 'data/champ_2000.dat' us qdidx:idx w dot, \
     'data/champ_2000_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2001 CHAMP Z residuals"
plot 'data/champ_2001.dat' us qdidx:idx w dot, \
     'data/champ_2001_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2002 CHAMP Z residuals"
plot 'data/champ_2002.dat' us qdidx:idx w dot, \
     'data/champ_2002_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2003 CHAMP Z residuals"
plot 'data/champ_2003.dat' us qdidx:idx w dot, \
     'data/champ_2003_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2004 CHAMP Z residuals"
plot 'data/champ_2004.dat' us qdidx:idx w dot, \
     'data/champ_2004_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2005 CHAMP Z residuals"
plot 'data/champ_2005.dat' us qdidx:idx w dot, \
     'data/champ_2005_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2006 CHAMP Z residuals"
plot 'data/champ_2006.dat' us qdidx:idx w dot, \
     'data/champ_2006_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2007 CHAMP Z residuals"
plot 'data/champ_2007.dat' us qdidx:idx w dot, \
     'data/champ_2007_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2008 CHAMP Z residuals"
plot 'data/champ_2008.dat' us qdidx:idx w dot, \
     'data/champ_2008_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "2009 CHAMP Z residuals"
plot 'data/champ_2009.dat' us qdidx:idx w dot, \
     'data/champ_2009_Z.dat' us 1:2 w li ls 3

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set xlabel "QD latitude (degrees)"
set title "2010 CHAMP Z residuals"
plot 'data/champ_2010.dat' us qdidx:idx w dot, \
     'data/champ_2010_Z.dat' us 1:2 w li ls 3
unset xlabel

unset multiplot
