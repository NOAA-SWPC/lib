#!/bin/sh
#
# Grid preproc output 'data.dat' file with season vs qdlat

datafile="data.dat"
outfile="data.grd"

qd_bin="5"
season_bin="20"

tmin="2000.0"
tmax="2006.0"

# maximum component value (nT)
compmax="40"

qdmax="70"

echo "time min = ${tmin}"
echo "time max = ${tmax}"

tidx="2"
qdidx="9"
compidx="19"

# Select qdlat in [-qdmax,qdmax], discard outliers
echo "selecting data from ${datafile}..."
tmp1=$(mktemp)
cat ${datafile} | \
    datasel -c ${qdidx} --min -${qdmax} --max ${qdmax} | \
    datasel -c ${tidx} --min ${tmin} --max ${tmax} | \
    datasel -c ${compidx} --neq 0.0 | \
    datasel -c 25 --eq 1 | \
    awk -v compidx=$compidx '{print $5,$9,$compidx}' > ${tmp1}

# Grid season,qdlat data
echo "gridding..."
cat ${tmp1} | gridgeo --xmin 0 --xmax 365 --ymin -${qdmax} --ymax ${qdmax} --xwidth ${season_bin} --ywidth ${qd_bin} -o ${outfile}

rm -f ${tmp1}

echo "output is ${outfile}"
