#!/bin/sh
#
# Generate ASCII/matlab coefficient files from binary

coefdir="/nfs/satmag_work/palken/corr/mfield/coef_new"
#outdir="/nfs/satmag_work/palken/corr/mfield/coef_ascii_new"
outdir="/nfs/satmag/SWARM/dmsp_ascii"
prog="/nfs/satmag_work/palken/corr/mfield/mfield_ascii"

rm -rf ${outdir}
mkdir ${outdir}

for file in $(ls ${coefdir}/*.dat); do
  echo "Processing ${file}..."
  outfile="${outdir}/$(basename ${file})"
  ${prog} -c ${file} -o ${outfile}
done

echo "Output directory is ${outdir}"
