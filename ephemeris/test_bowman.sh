#!/bin/sh
#
# Test Bowman ephemeris files

bowman_dir="/nfs/satmag_work/palken/dmsp_bowman"

prog="./test"

for sat in $(seq 15 18); do
  for year in $(seq 2009 2015); do
    file="${bowman_dir}/F${sat}_${year}_ephemeris.gz"

    if [[ -f "${file}" ]]; then
      echo "Testing ${file}..."
      ${prog} -b $file > /dev/null 2&>1
      rc=$?
      if [[ $rc != 0 ]]; then
        echo "FAILURE"
        exit
      fi
    fi
  done
done
