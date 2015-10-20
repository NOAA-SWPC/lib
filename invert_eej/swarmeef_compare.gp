#!/usr/bin/gnuplot

set term pngcairo enh col size 10in,8in

set style line 1 lt 1 lw 2
set style line 2 lt 2 lw 2
set style line 3 lt 3 lw 2

#set tmargin 0
#set bmargin 0
set lmargin 10
set rmargin 3

set xrange [2014.290341:2014.309302]

unset xtics

set out "swarmC_eef.png"
set multiplot layout 2,1 title "Swarm C EEF comparison"

set ylabel "EEF (mV/m)"
plot 'operC.dat' us 1:4 w li ti "OPER" ls 1, \
     'pproC.dat' us 1:4 w li ti "PPRO" ls 3

set xtics nomirror autofreq 0.003
set xlabel "time (years)"
set ylabel "percent error (%)"
plot 'oper_pproC.dat' us 1:(100*($4-$10)/$10) w li ti "% error" ls 2

unset multiplot

unset xtics

set out "swarmC_relerr.png"
set multiplot layout 2,1 title "Swarm C RelErr comparison"

set ylabel "RelErr"
plot 'operC.dat' us 1:5 w li ti "OPER" ls 1, \
     'pproC.dat' us 1:5 w li ti "PPRO" ls 3

set xtics nomirror autofreq 0.003
set xlabel "time (years)"
set ylabel "percent error (%)"
plot 'oper_pproC.dat' us 1:(100*($5-$11)/$11) w li ti "% error" ls 2

unset multiplot
