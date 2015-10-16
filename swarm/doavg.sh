#!/bin/sh
#
# Perform binning in QD latitude of vector residual data,
# output from dynamo.c
#
# Usage: doavg.sh [input_file]

input="dat"
if test -n "$1"; then
  input="$1"
fi

nbins="150"
ndata="500000"
qdmin="-55"
qdmax="55"

echo "input file = ${input}"

cat $input | datasel -c 7 --min ${qdmin} --max ${qdmax} | mbin -x 7 -y 8 -n ${nbins} -N ${ndata} > datX
cat $input | datasel -c 7 --min ${qdmin} --max ${qdmax} | mbin -x 7 -y 9 -n ${nbins} -N ${ndata} > datY
cat $input | datasel -c 7 --min ${qdmin} --max ${qdmax} | mbin -x 7 -y 10 -n ${nbins} -N ${ndata} > datZ

echo "Averaged outputs are: datX datY datZ"
