#!/bin/sh
#
# Plot temporal histogram of data coverage for model

prefix="datamap.dat"

# Number of days for each histogram bin
binsize="10"

dt=$(echo "scale=10; ${binsize} / 365.25" | bc)

echo "prefix = ${prefix}"
echo "bin size = ${binsize} [days]"

for sat in $(seq 0 2); do
  file="${prefix}.${sat}"
  if [ ! -f ${file} ]; then
    continue;
  fi

  echo "processing ${file}..."

  tmin=$(cat ${file} | datasel | awk '{print $1}' | head -1)
  tmax=$(cat ${file} | datasel | awk '{print $1}' | tail -1)

  nbins=$(echo "($tmax - $tmin) / $dt" | bc)

  outfile_scal="hist_scal.sat${sat}"
  outfile_vec="hist_vec.sat${sat}"
  outfile_euler="hist_euler.sat${sat}"

  echo "writing ${outfile_scal}..."
  cat ${file} | datasel -c 6 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_scal}

  echo "writing ${outfile_vec}..."
  cat ${file} | datasel -c 7 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_vec}

  echo "writing ${outfile_euler}..."
  cat ${file} | datasel -c 8 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_euler}
done
