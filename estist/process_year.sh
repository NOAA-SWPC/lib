#!/bin/sh
#
# Process 1 year of Est/Ist separation
# Usage: ./process_year.sh [year]

year="2000"
if test -n "$1"; then
  year="$1"
fi

outdir="$MYLIBHOME/estist/data"
prog="$MYLIBHOME/estist/main"

mkdir -p ${outdir}

# start/end times
t0=$(date -d "Jan 1 ${year} 00:00:00 UTC" +%s)
t1=$(date -d "Dec 31 ${year} 23:59:59 UTC" +%s)

outfile="${outdir}/Est_Ist_${year}.dat"
rm -f ${outfile}

echo "=== PROCESSING YEAR ${year} ==="
${prog} -a $t0 -b $t1 -o ${outfile}
