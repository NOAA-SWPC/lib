#!/home/palken/usr/bin/gnuplot

file1 = 'spectrum_Z2.s'
file2 = 'spectrum_F2.s'
file3 = 'diff.s'

nrow = 3
ncol = 1

load 'multi_default.cfg'

vbuffer = 1.2

load 'multi_defs.cfg'
load 'multi_png.cfg'
set out "spectrum.png"

load 'griddark.cfg'
load 'lines2.cfg'
set logscale y

set xlabel "spherical harmonic degree n"

set multiplot layout nrow,ncol

set title "Main field"
set ylabel "Power (nT^2)"
plot file1 us 1:2 w lp lw 4 ti "Model B", \
     file2 us 1:2 w lp lw 4 ti "Model C", \
     file3 us 1:2 w lp lw 4 ti "Difference"

load 'incrow.cfg'

set title "Secular variation"
set ylabel "Power (nT/year)^2"
plot [0:15] file1 us 1:($1 <= 15 ? $3 : 1/0) w lp lw 4 ti "Model B", \
            file2 us 1:($1 <= 15 ? $3 : 1/0) w lp lw 4 ti "Model C", \
            file3 us 1:($1 <= 15 ? $3 : 1/0) w lp lw 4 ti "Difference"

load 'incrow.cfg'

set title "Secular acceleration"
set ylabel "Power (nT/year^2)^2"
plot [0:8] file1 us 1:($1 <= 8 ? $4 : 1/0) w lp lw 4 ti "Model B", \
           file2 us 1:($1 <= 8 ? $4 : 1/0) w lp lw 4 ti "Model C", \
           file3 us 1:($1 <= 8 ? $4 : 1/0) w lp lw 4 ti "Difference"

unset multiplot
