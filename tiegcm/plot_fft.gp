#!/home/palken/usr/bin/gnuplot

set term pngcairo enh col size 1000,600

idx = 5
file = 'fft_data.txt'

set xlabel "days"
set ylabel "nT"

load 'grid.cfg'
load 'lines2.cfg'

set title "Windowed time series"

set out "fft1.png"
plot file us (column(0)/24):1 index idx w li lw 4 ti "q(2,0)"

set out "fft2.png"
plot file us (column(0)/24):1 index idx w li lw 4 ti "q(2,0)", \
     file us (column(0)/24):2 index idx w li lw 4 lt 3 ti "Hamming q(2,0)"

set xlabel "Period (days)"
set ylabel "Power (nT^2)"
set title "Windowed power spectrum"
set out "fft3.png"
plot file us 3:4 index idx w lp lw 4 lt 4 ti "q(2,0)"
