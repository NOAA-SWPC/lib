#!/bin/sh
#
# Make CHAMP index files

cdfdir="/data/CHAMP/Stage1_CHAOS"

start_year="2000"
end_year="2010"
start_month="01"
end_month="12"

# Make specialized index file
tmpfile=$(mktemp)
outfile="champ.idx"

for year in $(seq ${start_year} ${end_year}); do
  for month in $(seq -w ${start_month} ${end_month}); do
    find ${cdfdir} -name "CH-ME-3-MAG+${year}-${month}-*.cdf" -print >> ${tmpfile}
  done
done

echo "writing ${outfile}..."
cat ${tmpfile} | sort -g > ${outfile}
rm -f ${tmpfile}
