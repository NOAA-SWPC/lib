#!/usr/bin/gnuplot

set term pdfcairo enh col dashed

set out "extfield.pdf"

set style line 1 lt 1 lw 3 pt 3 lc rgb "red"
set style line 2 lt 1 lw 3 pt 3 lc rgb "green"
set style line 3 lt 1 lw 3 pt 3 lc rgb "blue"
set style line 4 lt 2 lw 3 pt 3 lc rgb "red"
set style line 5 lt 2 lw 3 pt 3 lc rgb "green"
set style line 6 lt 2 lw 3 pt 3 lc rgb "blue"
set style line 7 lt 1 lw 3 pt 3 lc rgb "black"
set style line 8 lt 1 lw 3 pt 3 lc rgb "purple"
set style line 9 lt 1 lw 3 pt 3 lc rgb "turquoise"

set key outside top horizontal font ",8" maxrows 1

set ylabel "nT"
set y2label "W/m^2/Hz"
set ytics nomirror
set y2tics

set y2range [-100:250]
set yrange [-30:50]

plot 'extfield.dat' us 1:2 w li ti "2014 k(1,0)" ls 1, \
     'extfield.dat' us 1:3 w li ti "2014 k(1,1)" ls 2, \
     'extfield.dat' us 1:4 w li ti "2014 k(1,-1)" ls 3, \
     'dat_ext.2013' us 1:2 w li ti "2013 k(1,0)" ls 4, \
     'dat_ext.2013' us 1:3 w li ti "2013 k(1,1)" ls 5, \
     'dat_ext.2013' us 1:4 w li ti "2013 k(1,-1)" ls 6, \
     'extfield.dat' us 1:7 w li ti "F107A (20 prev)" ls 7 axes x1y2, \
     'extfield.dat' us 1:5 w li ti "Est" ls 8, \
     'extfield.dat' us 1:6 w li ti "Ist" ls 9
