#!/usr/bin/env gnuplot

set term pngcairo enh col size 800,600
set out "current1d.png"

unset key

set cbrange [-80:80]
set cbtics 20

set cntrparam levels inc -80,10,80
#set cntrparam levels auto

lonA = -46.017235
lonB = -26.378850
lonC = -44.581345

ymin = -20
ymax = 20

idx=3

file = 'chi.txt'
tfile = 'test.dat'
cfile = 'cont.dat'

set table tfile
splot file us 1:2:3 index idx
unset table

set contour base
unset surface
set table cfile
splot file us 1:2:3 index idx
unset table
unset contour
set surface

load 'moreland.pal'
load 'xlonon.cfg'
load 'ylaton.cfg'
set cblabel "kA"

set arrow 1 from lonA,ymin to lonA,ymax front nohead lt 6 lw 2
#set arrow 2 from lonC,ymin to lonC,ymax front nohead lt 6 lw 2
#set arrow 3 from lonB,ymin to lonB,ymax front nohead lt 6 lw 2

#set label 1 "A,C" at lonA,ymax+1 center font ",8"
set label 1 "A" at lonA,ymax+1 center font ",8"
#set label 2 "B" at lonB,ymax+1 center font ",8"

set yrange [-20:20]
set ytics -20,5,20

#plot tfile w image, cfile w li lt -1 lw 1.5
plot tfile w image
