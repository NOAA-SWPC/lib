#!/usr/bin/gnuplot

set term pngcairo enh col size 10in,10in font ",18"
set out "dH.png"

set yrange [-0.8:0.8]

set xzeroaxis
set yzeroaxis

set key inside right bottom

set style line 1 lt 1 lw 2 lc rgb "black"

set title "Swarm EEF vs SAM-MBO {/Symbol \104}H measurements"
set xlabel "{/Symbol \104}H (nT)"
set ylabel "E_{SWARM} (mV/m)"
plot 'mboA' us 1:3 pt 5 lc rgb "red" ps 1.0 ti "Swarm A", \
     'mboB' us 1:3 pt 5 lc rgb "green" ps 1.0 ti "Swarm B", \
     'mboC' us 1:3 pt 5 lc rgb "blue" ps 1.0 ti "Swarm C"
