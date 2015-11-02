#!/bin/sh
#
# Plot maps of NEC data differences computed with two different
# Euler angles (fixed and 30-day angles)

prefix1="res_fixed"
prefix2="res_30d"

# A
sat="0"

time_start="2013.9"
time_end="2014.7"
time_step="0.1"

outdir="/nfs/satmag_work/palken/corr/mfield/nec_maps"

rm -rf ${outdir}
mkdir ${outdir}

function do_sat
{
  sat=$1
  t0=$2
  t1=$3

  outfile="${outdir}/nec_${t0}_sat${sat}"

  # X component difference
  paste ${prefix1}.sat${sat} ${prefix2}.sat${sat} | datasel -c 1 --min ${t0} --max ${t1} | awk '{print $4,$5,$11-$24}' > ${outfile}_X

  # Z component difference
  paste ${prefix1}.sat${sat} ${prefix2}.sat${sat} | datasel -c 1 --min ${t0} --max ${t1} | awk '{print $4,$5,$13-$26}' > ${outfile}_Z

  echo "Running gmtplot.sh..."
  sh gmtplot.sh ${outfile}_X
  sh gmtplot.sh ${outfile}_Z

  echo "Converting to png..."
  convert -flatten ${outfile}_X.ps ${outfile}_X.png
  convert -flatten ${outfile}_Z.ps ${outfile}_Z.png
}

for t0 in $(seq ${time_start} ${time_step} ${time_end}); do
  t1=$(echo ${t0} + ${time_step} | bc)

  echo "Processing time window [$t0,$t1]"

  do_sat 0 $t0 $t1
  do_sat 1 $t0 $t1
  do_sat 2 $t0 $t1

done

echo "Output directory is ${outdir}"
