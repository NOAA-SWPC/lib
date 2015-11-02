#!/usr/bin/gnuplot

set term pdfcairo enh col size 11in,8in font ",18"

set out "euler3.pdf"

set xlabel "time (years)"
set ylabel "angle (arcseconds)"

set xrange [2014.1:*]

set style line 1 lw 4
set style line 2 lw 4
set style line 3 lw 4

set multiplot layout 3,1

set title "Swarm A"
plot 'euler.comb.0' us 1:(($2-$6)*3600) w li ls 1 ti "Alpha", \
     'euler.comb.0' us 1:(($3-$7)*3600) w li ls 2 ti "Beta", \
     'euler.comb.0' us 1:(($4-$8)*3600) w li ls 3 ti "Gamma"

set title "Swarm B"
plot 'euler.comb.1' us 1:(($2-$6)*3600) w li ls 1 ti "Alpha", \
     'euler.comb.1' us 1:(($3-$7)*3600) w li ls 2 ti "Beta", \
     'euler.comb.1' us 1:(($4-$8)*3600) w li ls 3 ti "Gamma"

set title "Swarm C"
plot 'euler.comb.2' us 1:(($2-$6)*3600) w li ls 1 ti "Alpha", \
     'euler.comb.2' us 1:(($3-$7)*3600) w li ls 2 ti "Beta", \
     'euler.comb.2' us 1:(($4-$8)*3600) w li ls 3 ti "Gamma"

unset multiplot
