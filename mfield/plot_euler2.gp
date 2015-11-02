#!/usr/bin/gnuplot

set term pngcairo enh col size 20in,12in font ",10"

set out "euler2.png"

set xlabel "time (years)"
set ylabel "angle (degrees)"

set style line 1 lt 1 lw 2
set style line 2 lt 2 lw 2
set style line 4 lt 4 lw 2
set style line 5 lt 5 lw 2

#set xrange [2014.266:*]
#set xrange [2014.1:*]

set y2range [0:20000]

set ytics nomirror
set y2tics
set y2label "Vector data points per 10 days"

set xtics 0.2

set key font ",6" top left outside horizontal

set multiplot layout 3,3 columnsfirst

set title "Swarm A"

#set yrange [11.785:11.815]
plot 'hist_vec.sat0' us 1:3 w boxes fill so axes x1y2 lt 3 ti "", \
     'euler_6am6pm_30d.0' us 1:2 w li ls 1 ti "30d", \
     'euler_6am6pm_10d.0' us 1:2 w li ls 5 ti "10d", \
     11.8019 w li ls 4 ti "fixed", \
     11.8042 w li ls 2 ti "Nils"

#set yrange [-76.198:-76.192]
plot 'hist_vec.sat0' us 1:3 w boxes fill so axes x1y2 lt 3 ti "", \
     'euler_6am6pm_30d.0' us 1:3 w li ls 1 ti "30d", \
     'euler_6am6pm_10d.0' us 1:3 w li ls 5 ti "10d", \
     -76.1936 w li ls 4 ti "fixed", \
     -76.1930 w li ls 2 ti "Nils"

#set yrange [-12.59:-12.55]
plot 'hist_vec.sat0' us 1:3 w boxes fill so axes x1y2 lt 3 ti "", \
     'euler_6am6pm_30d.0' us 1:4 w li ls 1 ti "30d", \
     'euler_6am6pm_10d.0' us 1:4 w li ls 5 ti "10d", \
     -12.5692 w li ls 4 ti "fixed", \
     -12.5736 w li ls 2 ti "Nils"

set title "Swarm B"

#set yrange [-8.91:-8.87]
plot 'hist_vec.sat1' us 1:3 w boxes fill so axes x1y2 lt 3 ti "", \
     'euler_6am6pm_30d.1' us 1:2 w li ls 1 ti "30d", \
     'euler_6am6pm_10d.1' us 1:2 w li ls 5 ti "10d", \
     -8.8854 w li ls 4 ti "fixed", \
     -8.8847 w li ls 2 ti "Nils"

#set yrange [-76.431:-76.4265]
plot 'hist_vec.sat1' us 1:3 w boxes fill so axes x1y2 lt 3 ti "", \
     'euler_6am6pm_30d.1' us 1:3 w li ls 1 ti "30d", \
     'euler_6am6pm_10d.1' us 1:3 w li ls 5 ti "10d", \
     -76.4272 w li ls 4 ti "fixed", \
     -76.4266 w li ls 2 ti "Nils"

#set yrange [9.145:9.095]
plot 'hist_vec.sat1' us 1:3 w boxes fill so axes x1y2 lt 3 ti "", \
     'euler_6am6pm_30d.1' us 1:4 w li ls 1 ti "30d", \
     'euler_6am6pm_10d.1' us 1:4 w li ls 5 ti "10d", \
     9.1165 w li ls 4 ti "fixed", \
     9.1141 w li ls 2 ti "Nils"

set title "Swarm C"

#set yrange [1.795:1.83]
plot 'hist_vec.sat2' us 1:3 w boxes fill so axes x1y2 lt 3 ti "", \
     'euler_6am6pm_30d.2' us 1:2 w li ls 1 ti "30d", \
     'euler_6am6pm_10d.2' us 1:2 w li ls 5 ti "10d", \
     1.8120 w li ls 4 ti "fixed", \
     1.8130 w li ls 2 ti "Nils"

#set yrange [-76.833:-76.826]
plot 'hist_vec.sat2' us 1:3 w boxes fill so axes x1y2 lt 3 ti "", \
     'euler_6am6pm_30d.2' us 1:3 w li ls 1 ti "30d", \
     'euler_6am6pm_10d.2' us 1:3 w li ls 5 ti "10d", \
     -76.8273 w li ls 4 ti "fixed", \
     -76.8269 w li ls 2 ti "Nils"

#set yrange [-1.995:-1.96]
plot 'hist_vec.sat2' us 1:3 w boxes fill so axes x1y2 lt 3 ti "", \
     'euler_6am6pm_30d.2' us 1:4 w li ls 1 ti "30d", \
     'euler_6am6pm_10d.2' us 1:4 w li ls 5 ti "10d", \
     -1.9760 w li ls 4 ti "fixed", \
     -1.9789 w li ls 2 ti "Nils"

unset multiplot
