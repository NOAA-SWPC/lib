#!/bin/sh
#
# Plot temporal histogram of data coverage for model

prefix="output/datamap"

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
  outfile_vec_X="hist_vec_X.sat${sat}"
  outfile_vec_Y="hist_vec_Y.sat${sat}"
  outfile_vec_Z="hist_vec_Z.sat${sat}"
  outfile_euler="hist_euler.sat${sat}"

  echo "writing ${outfile_scal}..."
  cat ${file} | datasel -c 6 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_scal}

  echo "writing ${outfile_vec_X}..."
  cat ${file} | datasel -c 7 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_vec_X}

  echo "writing ${outfile_vec_Y}..."
  cat ${file} | datasel -c 8 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_vec_Y}

  echo "writing ${outfile_vec_Z}..."
  cat ${file} | datasel -c 9 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_vec_Z}

  echo "writing ${outfile_euler}..."
  cat ${file} | datasel -c 10 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_euler}
done
