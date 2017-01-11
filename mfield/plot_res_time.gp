#!/home/palken/usr/bin/gnuplot
#
# Plot time series of residuals from MF model
# for a single satellite. Also plotted is
# residuals vs QD latitude

file = 'res.dat.sat0'

nrow = 4
ncol = 2

load 'multi_default.cfg'

plotwidth = 8.0
b = 0.1
hbuffer = 1.0
r = -0.3

load 'multi_defs.cfg'
load 'multi_png.cfg'
set out "res_time.png"

set yrange [-30:30]
unset key

set multiplot layout nrow,ncol

set format x ""

set ylabel "X residual (nT)"
plot file us 1:($17 == 1 ? $9 : 1/0)

load 'incrow.cfg'

set ylabel "Y residual (nT)"
plot file us 1:($18 == 1 ? $10 : 1/0)

load 'incrow.cfg'

set yrange [-20:20]

set ylabel "Z residual (nT)"
plot file us 1:((abs($7) <= 55) & ($19 == 1) ? $11 : 1/0)

load 'incrow.cfg'

set xlabel "time (years)"
set format x "%g"

set yrange [*:*]
set ylabel "F residual (nT)"
plot file us 1:($16 == 1 ? $8 : 1/0)

load 'inccolumn.cfg'

load 'xlaton.cfg'
set format x ""
unset xlabel
set yrange [-30:30]

set ylabel "X residual (nT)"
plot file us 7:($17 == 1 ? $9 : 1/0)

load 'incrow.cfg'

set ylabel "Y residual (nT)"
plot file us 7:($18 == 1 ? $10 : 1/0)

load 'incrow.cfg'

set ylabel "Z residual (nT)"
plot file us 7:($19 == 1 ? $11 : 1/0)

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"
set format x "%g"

set yrange [*:*]
set ylabel "F residual (nT)"
plot file us 7:($16 == 1 ? $8 : 1/0)

unset multiplot
