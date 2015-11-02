#!/usr/bin/gnuplot

set term pngcairo enh col

set xlabel "QD latitude (degrees)"
set ylabel "residual (nT)"

set yrange [-20:20]
unset key

set out "Xscat.png"
set title "Swarm C X residuals"
plot [-55:55] 'res_30d_ext.sat2' us 6:8 w dot lt 3

set out "Yscat.png"
set title "Swarm C Y residuals"
plot [-55:55] 'res_30d_ext.sat2' us 6:9 w dot lt 3

set out "Zscat.png"
set title "Swarm C Z residuals"
plot [-55:55] 'res_30d_ext.sat2' us 6:10 w dot lt 3

set out "Fscat.png"
set title "Swarm C F residuals"
plot [-90:90] 'res_30d_ext.sat2' us 6:7 w dot lt 3
