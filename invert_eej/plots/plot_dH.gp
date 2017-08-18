#!/home/palken/usr/bin/gnuplot

nrow = 1
ncol = 2

load 'multi_default.cfg'

plotheight = 3.0
hbuffer = 1.0
r = -0.5
b = 0.3

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "plot_dH.png"

load 'xyborder.cfg'
load 'grid.cfg'

set xrange [-60:100]
set ylabel "peak current density (mA/m)"
set xlabel "{/Symbol \104}H (nT)"

set multiplot layout nrow,ncol

set key inside bottom right

set title "EEJ current vs JIC-PIU {/Symbol \104}H measurements"
plot 'jic_piu_A.dat' us 2:4 pt 5 lc rgb "red" ps 1.0 ti "Swarm A", \
     'jic_piu_B.dat' us 2:4 pt 5 lc rgb "green" ps 1.0 ti "Swarm B", \
     'jic_piu_C.dat' us 2:4 pt 5 lc rgb "blue" ps 1.0 ti "Swarm C"

unset key

load 'inccolumn.cfg'

set title "EEJ current vs SAM-MBO {/Symbol \104}H measurements"
plot 'sam_mbo_A.dat' us 2:4 pt 5 lc rgb "red" ps 1.0 ti "Swarm A", \
     'sam_mbo_B.dat' us 2:4 pt 5 lc rgb "green" ps 1.0 ti "Swarm B", \
     'sam_mbo_C.dat' us 2:4 pt 5 lc rgb "blue" ps 1.0 ti "Swarm C"

unset multiplot
