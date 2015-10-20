#!/bin/sh
# Make CHAMP index files for EEF processing

datadir="/nfs/satmag/CHAMP/Stage1"

for year in $(seq -w 00 10); do
  idxfile="champ${year}.idx"
  echo "Generating ${idxfile}..."
  find ${datadir} -name "CH-*MAG+20${year}-*.cdf" | sort -g > ${idxfile}
done
