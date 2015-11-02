#!/bin/sh
#
# Create index files suitable for modeling with 'mfield'
#
# Usage: domodel.sh [start_month] [start_year] [end_month] [end_year]
#
# Output files will be 'res_FXX.idx' for each satellite

start_month="1"
start_year="2011"
end_month="12"
end_year="2013"

if test -n "$1"; then
  start_month="$1"
fi
if test -n "$2"; then
  start_year="$2"
fi
if test -n "$3"; then
  end_month="$3"
fi
if test -n "$4"; then
  end_year="$4"
fi

suffix="BOW_FINAL"
dmsp_sats="F15 F16 F17 F18"

dmsp_dir="/nfs/satmag/DMSP/RES_CDF_${suffix}"

echo "start year = ${start_year}"
echo "end year = ${end_year}"

for s in ${dmsp_sats}; do
  # convert to lower case
  slow=$(echo ${s} | tr '[:upper:]' '[:lower:]')

  # initialize index file
  resfile="res_${slow}.idx"
  rm -f $resfile

  # construct index file
  echo "creating index ${resfile}"

  # add start_year starting at start_month
  for month in `seq -w ${start_month} 12`; do
    find ${dmsp_dir}/${start_year}/*${s}*${start_year}${month}*.cdf >> ${resfile}
  done

  # add middle years
  syp1=$((${start_year} + 1))
  eym1=$((${end_year} - 1))
  for year in `seq ${syp1} ${eym1}`; do
    find ${dmsp_dir}/${year}/*${s}*.cdf >> ${resfile}
  done

  # add end_year ending with end_month
  for month in `seq -f "%02g" 1 ${end_month}`; do
    find ${dmsp_dir}/${end_year}/*${s}*${end_year}${month}*.cdf >> ${resfile}
  done
done
