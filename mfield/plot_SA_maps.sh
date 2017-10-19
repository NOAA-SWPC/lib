#!/bin/sh

set term pngcairo enh col size 1000,1000

mapprog="$DATAHOME/palken/msynth-1.0/src/print_map"
plotprog="$DATAHOME/palken/msynth-1.0/src/plots/plotmap2.py"

coefdir="coef_F17"
outfile="F17.mp4"
title="DMSP F-17"

plot_args="-c "uT/yr^2" --cbmin -1.0 --cbmax 1.0 --cbtics 5 --cblev 1000"
#plot_args="-c "uT/yr^2""

# maximum SH degree for SA maps
nmax="6"

idx=1
for f in $(ls ${coefdir}/coef*.txt); do
  bname=$(basename $f ".txt")
  istr=$(seq -f "%02g" $idx $idx)
  outfile="${coefdir}/map.${istr}.png"
  epoch=$(echo $f | sed -r 's/.*\.([0-9]*\.[0-9]*)\..*/\1/g')
  echo "generating SA map for $f..."
  tmpfile=$(mktemp)
  ${mapprog} -c $f -n $nmax -o $tmpfile
  python ${plotprog} -i $tmpfile -o ${outfile} -t "${title}, epoch ${epoch}" ${plot_args}
  rm -f $tmpfile
  idx=$((idx+1))
done

#ffmpeg -framerate 4 -i map.%02d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf fps=25 ${outfile}
