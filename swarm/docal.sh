#!/bin/sh
#
# Do scalar calibration on Swarm data

datadir="/nfs/satmag/SWARM/TEST_EULER/DTU"
outdir="/nfs/satmag/SWARM/TEST_EULER/DTU_CAL"

param_file="params.dat"

prog="./testcal"

rm -rf $outdir
mkdir $outdir

rm -f ${param_file}
echo "# Field 1: time (years)" >> ${param_file}
echo "# Field 2: SX" >> ${param_file}
echo "# Field 3: SY" >> ${param_file}
echo "# Field 4: SZ" >> ${param_file}
echo "# Field 5: OX" >> ${param_file}
echo "# Field 6: OY" >> ${param_file}
echo "# Field 7: OZ" >> ${param_file}
echo "# Field 8: AXY" >> ${param_file}
echo "# Field 9: AXZ" >> ${param_file}
echo "# Field 10: AYZ" >> ${param_file}

for file in $(ls $datadir/*.cdf); do
  bname=$(basename ${file} .cdf)
  outfile="${outdir}/${bname}_CAL.cdf"

  ${prog} -i ${file} -o ${outfile} >> ${param_file}
done
