#!/bin/sh
#
# Make LT/QDlat grid of scalar data

datafile="data.dat"
outfile="data.grd"

# Component to grid (column of data file)
compidx="18"

min_lt="6"
max_lt="24"

lt_bin="1.5"
qd_bin="1.5"

tmin="2000.0"
tmax="2026.0"

# maximum component value (nT)
compmax="1000"

# maximum QD latitude (deg)
qdmax="70"

# min/max altitude (km)
altmin="0"
altmax="360"

echo "time min = ${tmin}"
echo "time max = ${tmax}"

echo "alt min = ${altmin}"
echo "alt max = ${altmax}"

echo "LT min = ${min_lt}"
echo "LT max = ${max_lt}"

tidx="2"
altidx="6"
qdidx="9"
gradidx="26"

# Select qdlat in [-qdmax,qdmax], discard outliers
echo "selecting data from ${datafile}..."
tmp1=$(mktemp)
cat ${datafile} | \
    datasel -c ${qdidx} --min -${qdmax} --max ${qdmax} | \
    datasel -c ${tidx} --min ${tmin} --max ${tmax} | \
    datasel -c ${altidx} --min ${altmin} --max ${altmax} | \
    datasel -c ${compidx} --min -${compmax} --max ${compmax} | \
    datasel -c ${gradidx} --eq 1 | \
    awk -v compidx=$compidx '{print $4,$9,$compidx}' > ${tmp1}

# Grid LT,qdlat data
echo "gridding..."
cat ${tmp1} | gridgeo --xmin ${min_lt} --xmax ${max_lt} --ymin -${qdmax} --ymax ${qdmax} --xwidth ${lt_bin} --ywidth ${qd_bin} -o ${outfile}

rm -f ${tmp1}

echo "output is ${outfile}"
