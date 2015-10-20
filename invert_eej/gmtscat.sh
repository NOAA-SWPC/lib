#!/bin/sh
#
# Make a GMT scatter plot of (x,y) data
# Usage: gmtscat.sh [xyfile]
#
# Where xyfile is formatted as: lon lat

infile="datamap.dat"
if test -n "$1"; then
  infile="$1"
fi

outfile="${infile}.ps"
rm -f $outfile

# psxy needs format: lat lon
xyfile=$(mktemp)
cat ${infile} | datasel | awk '{print $5,$4}' > ${xyfile}

# size of each scattered data point (circle)
pointsize="0.005i"

# Equatorial view
pscoast -X0.8i -Y6.5i -R-180/180/-90/90 -JW0/6.8 -B45g0f0/10g0f0 -N1 -Di -W2 -A10/1 -G180/180/180 -V -P -K >> $outfile
psxy $xyfile -: -R -JW -O -P -G255/0/0 -Sc${pointsize} -V -K >> $outfile

# North pole
pscoast -X0.4i -Y-2.6i -R-180/180/40/90 -JG0/90/2 -B45g0f0/10g0f0 -N1 -Di -W2 -A10/1 -G180/180/180 -V -P -O -K >> $outfile
psxy $xyfile -: -R -JG -O -P -G255/0/0 -Sc${pointsize} -V -K >> $outfile

# South pole
pscoast -X4.0i -Y0 -R-180/180/-90/-40 -JG0/-90/2 -B45g0f0/10g0f0 -N1 -Di -W2 -A10/1 -G180/180/180 -V -P -O -K >> $outfile
psxy $xyfile -: -R -JG -O -P -G255/0/0 -Sc${pointsize} -V >> $outfile

rm -f ${xyfile}

echo "Output is ${outfile}"
