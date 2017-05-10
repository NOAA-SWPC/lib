#!/bin/sh
#
# This script runs stage2 on DMSP data files to compute time
# shift and calibration parameters, and outputs a time series
# of each calibration parameter for a given satellite and year.
#
# Usage: stage2.sh [dmsp_sat] [year] [window_size]

dmsp_sat="F15"
if test -n "$1"; then
  dmsp_sat="$1"
fi

year="2008"
if test -n "$2"; then
  year="$2"
fi

# number of days in processing window (how many days
# to process at a time). Window is sliding
wsize="90"
if test -n "$3"; then
  wsize="$3"
fi

data_dir="$MYLIBHOME/dmsp/stage2_data"
mkdir -p ${data_dir}

# where to store model parameters
param_file="${data_dir}/params.${dmsp_sat}.${year}.w${wsize}.dat"

cdfdir="$DATAHOME/DMSP/MAG/Stage1"
prog="$MYLIBHOME/dmsp/stage2"

idxfile="$data_dir/${dmsp_sat}_${year}.idx"

ephfile="/nfs/satmag_work/palken/dmsp_bowman/${dmsp_sat}_${year}_ephemeris.gz"

function init_param_file
{
  pfile="$1"
  echo "parameter file = $pfile"
  rm -f $pfile
  echo "# Field 1: timestamp" >> $pfile
  echo "# Field 2: rms (nT)" >> $pfile
  echo "# Field 3: SX" >> $pfile
  echo "# Field 4: SY" >> $pfile
  echo "# Field 5: SZ" >> $pfile
  echo "# Field 6: OX (nT)" >> $pfile
  echo "# Field 7: OY (nT)" >> $pfile
  echo "# Field 8: OZ (nT)" >> $pfile
  echo "# Field 9: AXY (deg)" >> $pfile
  echo "# Field 10: AXZ (deg)" >> $pfile
  echo "# Field 11: AYZ (deg)" >> $pfile
}

echo "==== PROCESSING SATELLITE ${dmsp_sat} YEAR ${year} ===="

# Create index file of residuals
find ${cdfdir}/${year} -name "*${dmsp_sat}*.cdf" | sort -g > $idxfile

init_param_file "${param_file}"

i="1"
while [ $i -lt 13 ]; do
  t=$(echo "scale=8; $year + $i/12.0" | bc)
  tmin=$(echo "scale=8; $t - $wsize/365.25" | bc)
  tmax=$(echo "scale=8; $t + $wsize/365.25" | bc)

  # Calculate model parameters
  ${prog} -i $idxfile -p ${param_file} --t_min $tmin --t_max $tmax
  #echo "${prog} -i $idxfile -p ${param_file} --t_min $tmin --t_max $tmax"

  i=$[$i+1]
done

exit

for rstart in `seq 1 ${totfiles}`; do

  # Find number of data files to process based in wsize
  rend=`echo "${rstart} + ${wsize} - 1" | bc`
  if [ ${rend} -ge "${totfiles}" ]; then
    rend=${totfiles}
  fi

  tfile=`mktemp`

  # Pick out rows rstart to rend from idxfile
  cat ${idxfile} | awk "NR >= ${rstart} && NR <= ${rend}" > $tfile


  rm -f $tfile
done

rm -f ${idxfile}
