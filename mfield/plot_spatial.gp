#!/home/palken/usr/bin/gnuplot
#
# Plot time/latitude scatter plot of available data

nrow = 3
ncol = 3

load 'multi_default.cfg'

plotwidth = 6.0
b = 0.2
hbuffer = 1.0

load 'multi_defs.cfg'
load 'multi_png.cfg'

nskip = 3

set style fill transparent solid 1.0 noborder
set style circle radius 0.004

load 'lines2.cfg'
load 'ylaton.cfg'
#set xrange [2013.5:2016.5]
set xrange [2013.5:*]

satindices = "A B C"


do for [idx=0:2] {
  outfile = 'output/spatial_' . idx . '.png'
  set out outfile

  prefix = 'output/datamap' . idx
  outstr = sprintf('processing %s...', prefix)
  print outstr

  fileX = prefix . '_X.dat'
  fileY = prefix . '_Y.dat'
  fileZ = prefix . '_Z.dat'
  fileF = prefix . '_F.dat'
  fileDX_NS = prefix . '_DX_NS.dat'
  fileDY_NS = prefix . '_DY_NS.dat'
  fileDZ_NS = prefix . '_DZ_NS.dat'
  fileDX_EW = prefix . '_DX_EW.dat'
  fileDY_EW = prefix . '_DY_EW.dat'
  fileDZ_EW = prefix . '_DZ_EW.dat'

  load 'multi_reset.cfg'

  set multiplot layout nrow,ncol

  set format x ""
  unset key

  sat = word(satindices, idx + 1)
  set title 'Swarm '.sat

  set title "X component"
  plot fileX us 1:3 every nskip w dot ti "X"

  load 'incrow.cfg'

  set title "Y component"
  plot fileY us 1:3 every nskip w dot ti "Y"

  load 'incrow.cfg'

  set format x "%g"
  set xlabel "time"

  set title "Z component"
  plot fileZ us 1:3 every nskip w dot ti "Z"

  set format x ""
  unset xlabel

  load 'inccolumn.cfg'

  #plot fileF us 1:($3 > 0 ? $3 : 1/0) every nskip w dot ti "F"

  set title "N/S DX component"
  plot fileDX_NS us 1:3 every nskip w dot ti "N/S DX"

  load 'incrow.cfg'

  set title "N/S DY component"
  plot fileDY_NS us 1:3 every nskip w dot ti "N/S DY"

  load 'incrow.cfg'

  set format x "%g"
  set xlabel "time"

  set title "N/S DZ component"
  plot fileDZ_NS us 1:3 every nskip w dot ti "N/S DZ"

  set format x ""
  unset xlabel

  load 'inccolumn.cfg'

  set title "E/W DX component"
  plot fileDX_EW us 1:3 every nskip w dot ti "E/W DX"

  load 'incrow.cfg'

  set title "E/W DY component"
  plot fileDY_EW us 1:3 every nskip w dot ti "E/W DY"

  load 'incrow.cfg'

  set title "E/W DZ component"
  plot fileDZ_EW us 1:3 every nskip w dot ti "E/W DZ"

  unset multiplot

  outstr = sprintf('output is %s...', outfile)
  print outstr
}

