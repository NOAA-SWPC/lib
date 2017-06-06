#!/usr/bin/env gnuplot

set term pngcairo enh col size 1000,1000

coefdir = 'coef'
mapprog = '../msynth/print_map'

# maximum SH degree for SA maps
nmax = "6"

cmd = sprintf('ls %s/coef*.txt', coefdir)
list = system(cmd)

load 'xlonon.cfg'
load 'ylaton.cfg'
set pm3d map

set cblabel "uT/year^2"
set cbrange [-3:3]
set palette maxcolors 12

idx = 1
do for [file in list] {

  epoch = system(sprintf("echo %s | sed -e 's/[a-z,\/]*\.//' | sed -e 's/\.txt//'", file))

  outfile = sprintf('%s/map.%03d.txt', coefdir, idx)

  # Create map
  str = sprintf('Generating %s...', outfile)
  print str
  system(sprintf("%s -c %s -n %s -o %s", mapprog, file, nmax, outfile))

  pngfile = sprintf('%s/map.%03d.png', coefdir, idx)
  set out pngfile

  set title "Epoch ".epoch." degree ".nmax

  str = sprintf('Generating %s...', pngfile)
  print str
  splot outfile us 1:2:8

  idx = idx + 1
}
