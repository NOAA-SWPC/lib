#!/home/palken/usr/bin/gnuplot
#
# Plot time series (histogram) of model data distribution
# Use plothist.sh to generate data files

### n: change this parameter to equal the number of data sets to be plotted
nrow = 3
ncol = 1

load 'multi_default.cfg'
load 'multi_defs.cfg'
load 'multi_png.cfg'
set out "hist.png"

load 'lines2.cfg'

set ylabel "track scalar rms (nT)"
set xrange [2013.5:*]
set yrange [0:20000]
set ytics nomirror

set multiplot layout nrow,ncol

set xtics 0.5 format ""

set title "Temporal distribution of scalar data for main field modeling"
set ylabel "scalar data per 10 days"
plot 'hist_scal.sat0' us 1:3 w boxes fill so ti "Swarm A", \
     'hist_scal.sat1' us 1:3 w boxes fill so ti "Swarm B", \
     'hist_scal.sat2' us 1:3 w boxes fill so ti "Swarm C"

load 'incrow.cfg'

unset key

set title "Temporal distribution of vector data for main field modeling"
set ylabel "vector data per 10 days"
plot 'hist_vec.sat0' us 1:3 w boxes fill so, \
     'hist_vec.sat1' us 1:3 w boxes fill so, \
     'hist_vec.sat2' us 1:3 w boxes fill so

load 'incrow.cfg'

set xtics format "%g"

set title "Temporal distribution of vector data for Euler angles"
set ylabel "Euler angle data per 10 days"
plot 'hist_euler.sat0' us 1:3 w boxes fill so, \
     'hist_euler.sat1' us 1:3 w boxes fill so, \
     'hist_euler.sat2' us 1:3 w boxes fill so

unset multiplot
