#!/home/palken/usr/bin/gnuplot

set term pngcairo enh col size 1280,960 font ",22"

load 'xy2border.cfg'
load 'grid.cfg'
load 'lines2.cfg'

file = 'plp_high_night.dat'

set xrange [-60:60]

set yrange [-30:30]
set y2range [0:6*10**6]

set xlabel "QD latitude (degrees)"
set ylabel "residual (nT)"
set y2label "electron density (cm^{-3})"

set title "CHAMP, above 430 km altitude, 19-22 LT"
do for [idx=17:56] {
  set out sprintf('png/champ_max_night_%d.png', idx)
  plot file us 9:11 index idx w li lw 4 ti "X", file us 9:17 index idx w li axes x1y2 lw 4 ti "Ne"
}
