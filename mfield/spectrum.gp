#!/usr/bin/gnuplot

set term pngcairo enh col size 15in,15in font ",12"
set out "spectrum.png"

set style line 1 lw 2 lt 1
set style line 2 lw 2 lt 2
set style line 3 lw 2 lt 3

set logscale y

set xlabel "spherical harmonic degree"

set xrange [0:15]
set xtics autofreq 0,1,15

set multiplot layout 3,1

idx=2
set ylabel "power (nT^2)
set title "Main field spectrum"
plot 'wmm.s' us 1:idx w li ls 1 ti "WMM 2010", \
     'mfield.s' us 1:idx w li ls 2 ti "DMSP F16 2008-2010", \
     'pomme_spectra/model_76.spc' us 1:idx w li ls 3 ti "POMME"

idx=3
set ylabel "power (nT/year)^2
set title "Secular variation spectrum"
plot 'wmm.s' us 1:idx w li ls 1 ti "WMM 2010", \
     'mfield.s' us 1:idx w li ls 2 ti "DMSP F16 2008-2010", \
     'pomme_spectra/model_76.spc' us 1:idx w li ls 3 ti "POMME"

idx=4
set ylabel "power (nT/year^2)^2
set title "Secular acceleration spectrum"
plot 'wmm.s' us 1:idx w li ls 1 ti "WMM 2010", \
     'mfield.s' us 1:idx w li ls 2 ti "DMSP F16 2008-2010", \
     'pomme_spectra/model_76.spc' us 1:idx w li ls 3 ti "POMME"

unset multiplot
