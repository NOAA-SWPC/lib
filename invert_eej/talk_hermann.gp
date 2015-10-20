#!/usr/bin/gnuplot
# Figure for Hermann's EGU talk 2014

idx = 114

set style line 1 lw 3 lc rgb "red"
set style line 2 lw 3 lc rgb "blue"
set style line 3 lw 3 lc rgb "purple"
set style line 4 lw 3 lc rgb "green"

set xrange [-60:60]

set key inside bottom right
set term pngcairo enh col
set out "F1.png"
set xlabel "quasi-dipole latitude (degrees)"
set ylabel "scalar residual (nT)"
set title "Subtraction of main field model and Sq fit"
plot 'log/F2.dat' us 1:2 index idx w li ls 1 ti "scalar residual", \
     'log/F2.dat' us 1:3 index idx w li ls 4 ti "Sq model"

unset key
set term pngcairo enh col
set out "F2.png"
set xlabel "quasi-dipole latitude (degrees)"
set ylabel "scalar residual (nT)"
set title "Removal of Sq and external fields"
plot 'log/F2.dat' us 1:4 index idx w li ls 1

set xrange [-15:15]

set term pngcairo enh col
set out "EEJ.png"
set xlabel "quasi-dipole latitude (degrees)"
set ylabel "height-integrated eastward current density (A/m)"
set title "Inversion for EEJ current density"
plot 'log/EEJ.dat' us 1:2 index idx w li ls 2

set key inside top right

set term pngcairo enh col
set out "model.png"
set xlabel "quasi-dipole latitude (degrees)"
set ylabel "height-integrated eastward current density (A/m)"
set title "Modeling the EEJ current density"
plot 'log/model.dat' us 1:3 index idx w li ls 2 ti "Satellite", \
     'log/model.dat' us 1:2 index idx w li ls 3 ti "Model"

set term pngcairo enh col size 15in,12in
set out "model_compare.png"

unset ylabel
unset xlabel
unset title

set yrange [-0.1:0.15]

set multiplot layout 4,2

idx=56
set label 1 "RelErr = 0.63" at graph 0.01,0.9
plot 'log/model.dat' us 1:3 index idx w li ls 2 ti "Satellite", \
     'log/model.dat' us 1:2 index idx w li ls 3 ti "Model"
unset label 1

idx=102
set label 1 "RelErr = 0.61" at graph 0.01,0.9
plot 'log/model.dat' us 1:3 index idx w li ls 2 ti "Satellite", \
     'log/model.dat' us 1:2 index idx w li ls 3 ti "Model"
unset label 1

idx=112
set label 1 "RelErr = 0.55" at graph 0.01,0.9
plot 'log/model.dat' us 1:3 index idx w li ls 2 ti "Satellite", \
     'log/model.dat' us 1:2 index idx w li ls 3 ti "Model"
unset label 1

idx=132
set label 1 "RelErr = 0.84" at graph 0.01,0.9
plot 'log/model.dat' us 1:3 index idx w li ls 2 ti "Satellite", \
     'log/model.dat' us 1:2 index idx w li ls 3 ti "Model"
unset label 1

idx=142
set label 1 "RelErr = 0.43" at graph 0.01,0.9
plot 'log/model.dat' us 1:3 index idx w li ls 2 ti "Satellite", \
     'log/model.dat' us 1:2 index idx w li ls 3 ti "Model"
unset label 1

idx=6
set label 1 "RelErr = 0.57" at graph 0.01,0.9
plot 'log/model.dat' us 1:3 index idx w li ls 2 ti "Satellite", \
     'log/model.dat' us 1:2 index idx w li ls 3 ti "Model"
unset label 1

set xlabel "quasi-dipole latitude (degrees)"
idx=162
set label 1 "RelErr = 0.67" at graph 0.01,0.9
plot 'log/model.dat' us 1:3 index idx w li ls 2 ti "Satellite", \
     'log/model.dat' us 1:2 index idx w li ls 3 ti "Model"
unset label 1

idx=120
set label 1 "RelErr = 0.76" at graph 0.01,0.9
plot 'log/model.dat' us 1:3 index idx w li ls 2 ti "Satellite", \
     'log/model.dat' us 1:2 index idx w li ls 3 ti "Model"
unset label 1

unset multiplot
