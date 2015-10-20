#!/bin/sh
#
# Run EEF chain on CHAMP data

prog="./main"

# set for profiles only
flags="-p"
#flags=""

for year in $(seq -w 00 05); do
  log_dir="log_champ_${year}"
  idx_file="champ${year}.idx"

  rm -rf ${log_dir}
  mkdir ${log_dir}

  echo "Processing year ${year} (${idx_file}, ${log_dir})"
  screen -d -m -S champ${year} ${prog} -c ${idx_file} -l ${log_dir} ${flags}
done
