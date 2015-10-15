#!/bin/sh
# Print time series of CHAMP mean altitude

sh mkidx.sh

outfile="champ_alt.dat"
rm -f ${outfile}

prog="./print_alt"

idxfile=$(mktemp)
cat champ_20*.idx | sort -g > ${idxfile}

${prog} -i ${idxfile} > ${outfile}

rm -f ${idxfile}
echo "Output is ${outfile}"
