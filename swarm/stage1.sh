#!/bin/sh
#
# This script computes main, crustal, and external field
# values for Swarm data
#
# It is designed to be run continually, so as new data
# files are downloaded, the script checks if a file has
# already been processed before computing the fields
#
# Usage: ./stage1.sh [A|B|C] [year]

# Keep every nth sample
down_sample="1"

cdfdir="/data/SWARM/MAG/Unzipped_Data"

outdir="/data/SWARM/MAG/Stage1"
#outdir="/data/SWARM/MAG/Stage1_CHAOS"

# Use CHAOS external field model
#extra_flags="-c"
extra_flags=""

prog="/data/palken/lib/swarm/stage1"

sat=""
if test -n "$1"; then
  sat="$1"
fi

year=""
if test -n "$2"; then
  year="$2"
fi

#rm -rf ${outdir}
mkdir ${outdir}

for file in $(ls ${cdfdir}/*MAG${sat}*LR_1B_${year}*.cdf); do
  bname=$(basename "${file}" .cdf)
  outfile="${outdir}/${bname}_Stage1.cdf"

  # Check first if we already processed this file
  if [[ ! -f "${outfile}" ]]; then
    echo "Processing: ${file}"
    ${prog} -d ${down_sample} -i ${file} -o ${outfile} ${extra_flags}
  fi
done
