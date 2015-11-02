#!/bin/sh
#
# Make GMT lat/lon maps of the X,Y,Z,F residuals from a model run

file="res.sat0"

tmpfile=$(mktemp)

# select QD latitude < 55
cat $file | datasel -c 6 --min -55 --max 55 > $tmpfile

cat $tmpfile | awk '{print $4,$5,$8}' > dat
sh gmtplot.sh dat
convert -flatten dat.ps X.png

cat $tmpfile | awk '{print $4,$5,$9}' > dat
sh gmtplot.sh dat
convert -flatten dat.ps Y.png

cat $tmpfile | awk '{print $4,$5,$10}' > dat
sh gmtplot.sh dat
convert -flatten dat.ps Z.png

cat $file | datasel | awk '{print $4,$5,$7}' > dat
sh gmtplot.sh dat
convert -flatten dat.ps F.png

rm -f $tmpfile
