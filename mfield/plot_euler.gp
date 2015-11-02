#!/usr/bin/gnuplot

### n: change this parameter to equal the number of data sets to be plotted
nrow = 5
ncol = 3

# Height of each plot in pixels
plotheight = 150.0

# Width of each plot in pixels
plotwidth = 400.0

# Vertical buffer space between each plot in pixels
vbuffer = 20

# Horizontal buffer space between each plot in pixels
hbuffer = 80

# t: top margin in pixels
t = 75.0

# b: key height in pixels (bottom margin)
b = 200.0

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
set out "euler.png"

set ylabel "angle (degrees)"
set xrange [2013.9:2014.9]

set style line 1 lt 1 lw 2
set style line 2 lt 2 lw 2
set style line 3 lc rgb "blue" lw 2
set style line 4 lt 4 lw 2
set style line 5 lt 5 lw 2
set style line 6 lc rgb "orange" lw 2
set style line 7 lc rgb "brown" lw 2
set style line 8 lc rgb "black" lw 2

set xtics 0.2 format "" nomirror

# Offsets to add to swarm A curves in arcseconds
offset_A_alpha = 200.0
offset_A_gamma = -200.0

preflight_A_alpha = 11.8294 - offset_A_alpha/3600
preflight_A_beta = -76.1904
preflight_A_gamma = -12.5941 - offset_A_gamma/3600
preflight_B_alpha = -8.9132
preflight_B_beta = -76.4241
preflight_B_gamma = 9.1525
preflight_C_alpha = 1.7931
preflight_C_beta = -76.8270
preflight_C_gamma = -1.9536


set multiplot layout (nrow+1),ncol columnsfirst

set lmargin at screen left(1,ncol,w,l,r,hbuffer)
set rmargin at screen right(1,ncol,w,l,r,hbuffer)

unset key

# First plot

currentplot = 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set title "Swarm A"

set ylabel "{/Symbol \104}{/Symbol \141} (arcseconds)"
set yrange [-20:180]
plot (11.8019 - preflight_A_alpha)*3600 w li ls 4 ti "Fixed", \
     (11.8042 - preflight_A_alpha)*3600 w li ls 2 ti "DTU", \
     'euler_6am6pm_10d.0' us 1:($2-preflight_A_alpha)*3600 w li ls 5 ti "10d", \
     'euler_6am6pm_30d.0' us 1:($2-preflight_A_alpha)*3600 w li ls 1 ti "30d", \
     'euler_6am6pm_60d.0' us 1:($2-preflight_A_alpha)*3600 w li ls 6 ti "60d"
set yrange [*:*]

unset title

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \142} (arcseconds)"
set yrange [-30:15]
plot (-76.1936 - preflight_A_beta)*3600 w li ls 4 ti "Fixed", \
     (-76.1930 - preflight_A_beta)*3600 w li ls 2 ti "DTU", \
     'euler_6am6pm_10d.0' us 1:($3-preflight_A_beta)*3600 w li ls 5 ti "10d", \
     'euler_6am6pm_30d.0' us 1:($3-preflight_A_beta)*3600 w li ls 1 ti "30d", \
     'euler_6am6pm_60d.0' us 1:($3-preflight_A_beta)*3600 w li ls 6 ti "60d"
set yrange [*:*]

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \147} (arcseconds)"
set yrange [-250:50]
plot (-12.5692 - preflight_A_gamma)*3600 w li ls 4 ti "Fixed", \
     (-12.5736 - preflight_A_gamma)*3600 w li ls 2 ti "DTU", \
     'euler_6am6pm_10d.0' us 1:($4-preflight_A_gamma)*3600 w li ls 5 ti "10d", \
     'euler_6am6pm_30d.0' us 1:($4-preflight_A_gamma)*3600 w li ls 1 ti "30d", \
     'euler_6am6pm_60d.0' us 1:($4-preflight_A_gamma)*3600 w li ls 6 ti "60d"
set yrange [*:*]

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "local time (hours)"
set yrange [0:24]
set ytics 4
set y2range [-5:5]
set label 1 "north" at graph 0.65,0.6
set label 2 "south" at graph 0.65,0.45
plot 'res_10d_qdfilt.sat0' us 1:2 pt 6 ps 0.5 lc rgb "black", \
     6 ls 8, \
     18 ls 8, \
     'res_10d_qdfilt.sat0' us 1:14 axes x1y2 w li ls 7
set yrange [*:*]
set ytics auto
set y2range [*:*]
unset label 1
unset label 2

## Last plot

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)
set xlabel "time (years)"
set xtics format "%g"

set ylabel "vector data / 10 days"
set yrange [0:14000]
plot 'hist_vec.sat0' us 1:3 w boxes fill so lt 3 ti ""
set yrange [*:*]

# Key plot
set tmargin at screen bot(nrow,nrow,h,t,b,vbuffer)
set bmargin at screen 0
#set key center center horizontal
set border 0
unset tics
unset xlabel
unset ylabel
set yrange [0:1]

plot 2 ls 4 t "Fixed", \
     2 ls 2 t "DTU/0302", \
     2 ls 5 t "10-day", \
     2 ls 1 t "30-day"

set title "Swarm B"

unset key
set yrange [*:*]
set border
set ytics
set xtics 0.2 format "" nomirror
set lmargin at screen left(2,ncol,w,l,r,hbuffer)
set rmargin at screen right(2,ncol,w,l,r,hbuffer)

currentplot = 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \141} (arcseconds)"
set yrange [-20:180]
plot (-8.8854 - preflight_B_alpha)*3600 w li ls 4 ti "Fixed", \
     (-8.8847 - preflight_B_alpha)*3600 w li ls 2 ti "DTU", \
     'euler_6am6pm_10d.1' us 1:($2-preflight_B_alpha)*3600 w li ls 5 ti "10d", \
     'euler_6am6pm_30d.1' us 1:($2-preflight_B_alpha)*3600 w li ls 1 ti "30d", \
     'euler_6am6pm_60d.1' us 1:($2-preflight_B_alpha)*3600 w li ls 6 ti "60d"
set yrange [*:*]

unset title

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \142} (arcseconds)"
set yrange [-30:15]
plot (-76.4272 - preflight_B_beta)*3600 w li ls 4 ti "Fixed", \
     (-76.4266 - preflight_B_beta)*3600 w li ls 2 ti "DTU", \
     'euler_6am6pm_10d.1' us 1:($3-preflight_B_beta)*3600 w li ls 5 ti "10d", \
     'euler_6am6pm_30d.1' us 1:($3-preflight_B_beta)*3600 w li ls 1 ti "30d", \
     'euler_6am6pm_60d.1' us 1:($3-preflight_B_beta)*3600 w li ls 6 ti "60d"
set yrange [*:*]

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \147} (arcseconds)"
set yrange [-250:50]
plot (9.1165 - preflight_B_gamma)*3600 w li ls 4 ti "Fixed", \
     (9.1141 - preflight_B_gamma)*3600 w li ls 2 ti "DTU", \
     'euler_6am6pm_10d.1' us 1:($4-preflight_B_gamma)*3600 w li ls 5 ti "10d", \
     'euler_6am6pm_30d.1' us 1:($4-preflight_B_gamma)*3600 w li ls 1 ti "30d", \
     'euler_6am6pm_60d.1' us 1:($4-preflight_B_gamma)*3600 w li ls 6 ti "60d"
set yrange [*:*]

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)
set ylabel "local time (hours)"
set yrange [0:24]
set ytics 4
set y2range [-5:5]
set label 1 "north" at graph 0.67,0.6
set label 2 "south" at graph 0.67,0.45
plot 'res_10d_qdfilt.sat1' us 1:2 pt 6 ps 0.5 lc rgb "black", \
     6 ls 8, \
     18 ls 8, \
     'res_10d_qdfilt.sat1' us 1:14 axes x1y2 w li ls 7
set yrange [*:*]
set ytics auto
set y2range [*:*]
unset label 1
unset label 2

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)
set xlabel "time (years)"
set xtics format "%g"
set ylabel "vector data / 10 days"
set yrange [0:14000]
plot 'hist_vec.sat1' us 1:3 w boxes fill so lt 3 ti ""
set yrange [*:*]

# Key plot
set tmargin at screen bot(nrow,nrow,h,t,b,vbuffer)
set bmargin at screen 0
set key center center horizontal
set border 0
unset tics
unset xlabel
unset ylabel
set yrange [0:1]

plot 2 ls 4 t "Fixed - preflight", \
     2 ls 2 t "DTU/0302 - preflight", \
     2 ls 5 t "10-day - preflight", \
     2 ls 1 t "30-day - preflight", \
     2 ls 6 t "60-day - preflight"

set title "Swarm C"

unset key
set yrange [*:*]
set border
set ytics
set xtics 0.2 format "" nomirror
set lmargin at screen left(3,ncol,w,l,r,hbuffer)
set rmargin at screen right(3,ncol,w,l,r,hbuffer)

currentplot = 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \141} (arcseconds)"
set yrange [-20:180]
plot (1.8120 - preflight_C_alpha)*3600 w li ls 4 ti "Fixed", \
     (1.8130 - preflight_C_alpha)*3600 w li ls 2 ti "DTU", \
     'euler_6am6pm_10d.2' us 1:($2-preflight_C_alpha)*3600 w li ls 5 ti "10d", \
     'euler_6am6pm_30d.2' us 1:($2-preflight_C_alpha)*3600 w li ls 1 ti "30d", \
     'euler_6am6pm_60d.2' us 1:($2-preflight_C_alpha)*3600 w li ls 6 ti "60d"
set yrange [*:*]

unset title

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \142} (arcseconds)"
set yrange [-30:15]
plot (-76.8273 - preflight_C_beta)*3600 w li ls 4 ti "Fixed", \
     (-76.8269 - preflight_C_beta)*3600 w li ls 2 ti "DTU", \
     'euler_6am6pm_10d.2' us 1:($3-preflight_C_beta)*3600 w li ls 5 ti "10d", \
     'euler_6am6pm_30d.2' us 1:($3-preflight_C_beta)*3600 w li ls 1 ti "30d", \
     'euler_6am6pm_60d.2' us 1:($3-preflight_C_beta)*3600 w li ls 6 ti "60d"
set yrange [*:*]

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)

set ylabel "{/Symbol \104}{/Symbol \147} (arcseconds)"
set yrange [-250:50]
plot (-1.9760 - preflight_C_gamma)*3600 w li ls 4 ti "Fixed", \
     (-1.9789 - preflight_C_gamma)*3600 w li ls 2 ti "DTU", \
     'euler_6am6pm_10d.2' us 1:($4-preflight_C_gamma)*3600 w li ls 5 ti "10d", \
     'euler_6am6pm_30d.2' us 1:($4-preflight_C_gamma)*3600 w li ls 1 ti "30d", \
     'euler_6am6pm_60d.2' us 1:($4-preflight_C_gamma)*3600 w li ls 6 ti "60d"
set yrange [*:*]

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)
set ylabel "local time (hours)"
set yrange [0:24]
set ytics 4
set y2range [-5:5]
set label 1 "north" at graph 0.65,0.6
set label 2 "south" at graph 0.65,0.45
plot 'res_10d_qdfilt.sat2' us 1:2 pt 6 ps 0.5 lc rgb "black", \
     6 ls 8, \
     18 ls 8, \
     'res_10d_qdfilt.sat2' us 1:14 axes x1y2 w li ls 7
set yrange [*:*]
set ytics auto
set y2range [*:*]
unset label 1
unset label 2

currentplot = currentplot + 1
set tmargin at screen top(currentplot,nrow,h,t,b,vbuffer)
set bmargin at screen bot(currentplot,nrow,h,t,b,vbuffer)
set xlabel "time (years)"
set xtics format "%g"
set ylabel "vector data / 10 days"
set yrange [0:14000]
plot 'hist_vec.sat2' us 1:3 w boxes fill so lt 3 ti ""
set yrange [*:*]

# Key plot
set tmargin at screen bot(nrow,nrow,h,t,b,vbuffer)
set bmargin at screen 0
#set key center center horizontal
set border 0
unset tics
unset xlabel
unset ylabel
set yrange [0:1]

plot 2 ls 4 t "Fixed", \
     2 ls 2 t "DTU/0302", \
     2 ls 5 t "10-day", \
     2 ls 1 t "30-day"

unset multiplot
