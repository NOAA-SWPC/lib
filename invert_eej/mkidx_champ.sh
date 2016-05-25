#!/bin/sh
# Make CHAMP index files for EEF processing

datadir="$DATAHOME/CHAMP/Stage1_CHAOS"

for year in $(seq -w 00 10); do
  idxfile="champ${year}.idx"
  echo "Generating ${idxfile}..."
  find ${datadir} -name "CH-*MAG+20${year}-*.cdf" | sort -g > ${idxfile}
done
