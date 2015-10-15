#!/bin/sh
#
# This script computes main,crustal,external field values for CHAMP data
#
# Usage: ./stage1.sh [year]

# Keep every nth sample
down_sample="1"

cdfdir="/nfs/satmag/CHAMP/L3"
#outdir="/nfs/satmag/CHAMP/Stage1"
outdir="/nfs/satmag/CHAMP/Stage1_CHAOS"

#cdfdir="/nfs/satmag/CHAMP/test"
#outdir="/nfs/satmag/CHAMP/Stage1_TEST"

prog="/nfs/satmag_work/palken/corr/champ/stage1"

# Use CHAOS external field model?
extra_flags="-c"

year=""
if test -n "$1"; then
  year="$1"
fi

#rm -rf ${outdir}
mkdir ${outdir}

if [ -z "${year}" ]; then
  yearlist=$(seq 2000 2010)
else
  yearlist="$year"
fi

for yr in ${yearlist}; do
  filelist=$(ls ${cdfdir}/${yr}/*.cdf)
#  filelist=$(ls ${cdfdir}/*.cdf)

  for file in ${filelist}; do
    bname=$(basename "${file}" .cdf)
    outfile="${outdir}/${bname}_Stage1.cdf"

    # Check first if we already processed this file
    if [[ ! -f "${outfile}" ]]; then
      echo "Processing: ${file}"
      ${prog} -d ${down_sample} -i ${file} -o ${outfile} ${extra_flags}
    fi
  done
done
