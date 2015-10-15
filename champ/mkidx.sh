#!/bin/sh
#
# Make CHAMP index files
#
# Usage: mkidx.sh [start_year] [end_year] [start_doy] [end_doy] [start_doy2] [end_doy2]

cdfdir="/nfs/satmag/CHAMP/Stage1"

start_year="2000"
end_year="2010"
start_doy="1"
end_doy="366"
start_doy2=""
end_doy2=""

if test -n "$1"; then start_year="$1"; fi
if test -n "$2"; then end_year="$2"; fi
if test -n "$3"; then start_doy=$(seq -f "%02g" $3 $3); fi
if test -n "$4"; then end_doy=$(seq -f "%02g" $4 $4); fi
if test -n "$5"; then start_doy2=$(seq -f "%02g" $5 $5); fi
if test -n "$6"; then end_doy2=$(seq -f "%02g" $6 $6); fi

doy_prog="/nfs/satmag_work/palken/corr/common/doy2md"

# Make specialized index file
tmpfile=$(mktemp)
outfile="champ.idx"

for year in $(seq ${start_year} ${end_year}); do

  for doy in $(seq ${start_doy} ${end_doy}); do
    month=$(${doy_prog} ${year} ${doy} | awk '{print $1}')
    day=$(${doy_prog} ${year} ${doy} | awk '{print $2}')
    find ${cdfdir} -name "CH-ME-3-MAG+${year}-${month}-${day}*.cdf" -print >> ${tmpfile}
  done

  # add additional days if requested
  if [ -n "${start_doy2}" ] && [ -n "${end_doy2}" ]; then
    for doy in $(seq ${start_doy2} ${end_doy2}); do
      month=$(${doy_prog} ${year} ${doy} | awk '{print $1}')
      day=$(${doy_prog} ${year} ${doy} | awk '{print $2}')
      find ${cdfdir} -name "CH-ME-3-MAG+${year}-${month}-${day}*.cdf" -print >> ${tmpfile}
    done
  fi
done

echo "writing ${outfile}..."
cat ${tmpfile} | sort -g > ${outfile}
rm -f ${tmpfile}

# Make yearly index files
for year in $(seq 2000 2010); do
  outfile="champ_${year}.idx"
  echo "writing ${outfile}..."
  find ${cdfdir} -name "CH-ME-3-MAG+${year}*.cdf" -print | sort -g > ${outfile}
done
