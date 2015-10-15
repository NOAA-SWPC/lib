#!/bin/sh
#
# Plot CHAMP residuals to search for gravity and diamagnetic signatures

# Where to store data files
outdir="/nfs/satmag_work/palken/corr/champ/data"

# Local time range
lt_min="20"
lt_max="21"

# QD latitude range
qd_min="-50"
qd_max="50"

# Altitude range
alt_min="0"
alt_max="1000"

# Max KP
kp_max="2"

# Downsample factor
downsample="1"

start_doy="1"
end_doy="365"

#start_doy2="11"
#end_doy2="12"
start_doy2=""
end_doy2=""

rm -rf ${outdir}
mkdir -p ${outdir}

for year in $(seq 2000 2010); do
  # Make champ.idx file with selected year and months
  sh mkidx.sh ${year} ${year} ${start_doy} ${end_doy} ${start_doy2} ${end_doy2}

  outfile="${outdir}/champ_${year}.dat"
  echo "Generating ${outfile}..."
  ./print -i champ.idx --lt_min ${lt_min} --lt_max ${lt_max} --alt_min ${alt_min} --alt_max ${alt_max} --qd_min ${qd_min} --qd_max ${qd_max} --kp_max ${kp_max} -d ${downsample} > ${outfile}

  avgfile="${outdir}/champ_${year}_X.dat"
  echo "Generating ${avgfile}..."
  cat ${outfile} | mbin -x 10 -y 17 -n 40 -N 2000000 > ${avgfile}

  avgfile="${outdir}/champ_${year}_Y.dat"
  echo "Generating ${avgfile}..."
  cat ${outfile} | mbin -x 10 -y 18 -n 40 -N 2000000 > ${avgfile}

  avgfile="${outdir}/champ_${year}_Z.dat"
  echo "Generating ${avgfile}..."
  cat ${outfile} | mbin -x 10 -y 19 -n 40 -N 2000000 > ${avgfile}
done
