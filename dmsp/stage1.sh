#!/bin/sh
#
# This script computes main field values for DMSP data
#
# Usage: stage1.sh [FXX] [year]

# satellite to process
dmsp_sat="F15"
if test -n "$1"; then
  dmsp_sat="$1"
fi

# year to process
year="2009"
if test -n "$2"; then
  year="$2"
fi

# Keep every nth sample
down_sample="1"

indir="$DATAHOME/DMSP/MAG/Stage0"
outdir="$DATAHOME/DMSP/MAG/Stage1"

prog="$MYLIBHOME/dmsp/stage1"

# Extra flags (such as use CHAOS)
extra_flags="-c"

echo "==== PROCESSING YEAR ${year} ===="

mkdir -p $outdir/$year

for file in $(ls ${indir}/${year}/DMSP_${dmsp_sat}_*.cdf); do
  bname=$(basename "${outdir}/${year}/${file}" _Stage0.cdf)
  outfile="${outdir}/${year}/${bname}_Stage1.cdf"

  # Check first if we already processed this file
  if [[ ! -f "${outfile}" ]]; then
    echo "Processing: ${file}"
    ${prog} -d ${down_sample} -i ${file} -o ${outfile} ${extra_flags}
  fi
done
