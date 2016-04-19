#!/usr/bin/env gnuplot
#
# Plot magnetic field response for first few PCs

nrow = 3
ncol = 3

load 'multi_default.cfg'

plotwidth = 5.0
plotheight = 2.0
b = 0.3
hbuffer = 1.5

load 'multi_defs.cfg'
load 'multi_png.cfg'

set output "B.png"

file = 'pc_time.txt'

comps = "B_x B_y B_z"

unset key

set xrange [-180:180]
set xtics -180,45,180

set yrange [-75:75]
set ytics -75,15,75

set ylabel "latitude (degrees)"
set cblabel "dimensionless"

set pm3d map
load 'moreland.pal'

set multiplot layout nrow,ncol

pc1_idx = 4
pc2_idx = 8
pc3_idx = 12

do for [idx=1:3] {

  set format x ""
  unset xlabel

  if (idx == 1) {
    set cbrange [-10:10]
  } else {
    if (idx == 2) {
      set cbrange [-4:4]
    } else {
      set cbrange [-10:10]
    }
  }

  set title word(comps,idx)
  splot file us 1:2:pc1_idx
  unset title

  load 'incrow.cfg'

  splot file us 1:2:pc2_idx

  load 'incrow.cfg'

  set xlabel "longitude (degrees)"
  set format x "%g"

  splot file us 1:2:pc3_idx

  load 'inccolumn.cfg'

  pc1_idx = pc1_idx + 1
  pc2_idx = pc2_idx + 1
  pc3_idx = pc3_idx + 1
}

unset multiplot
