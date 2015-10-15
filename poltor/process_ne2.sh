#!/bin/sh
#
# Grid preproc output 'data.dat' file with LT vs qdlat

datafile="data.dat"
outfile="data.grd"

lt_bin="1.5"
qd_bin="2"

tmin="2000.0"
tmax="2026.0"

# maximum component value (nT)
compmax="40"

qdmax="70"

echo "time min = ${tmin}"
echo "time max = ${tmax}"

tidx="2"
ltidx="4"
qdidx="9"
compidx="19"

# Select qdlat in [-qdmax,qdmax], discard outliers
echo "selecting data from ${datafile}..."
tmp1=$(mktemp)
cat ${datafile} | \
    datasel -c ${qdidx} --min -${qdmax} --max ${qdmax} | \
    datasel -c ${tidx} --min ${tmin} --max ${tmax} | \
    datasel -c ${ltidx} --min 6 --max 24 | \
    datasel -c ${compidx} --neq 0.0 | \
    awk -v compidx=$compidx '{print $4,$9,$compidx}' > ${tmp1}

# Grid season,qdlat data
echo "gridding..."
cat ${tmp1} | gridgeo --xmin 6 --xmax 24 --ymin -${qdmax} --ymax ${qdmax} --xwidth ${lt_bin} --ywidth ${qd_bin} -o ${outfile}

rm -f ${tmp1}

echo "output is ${outfile}"
