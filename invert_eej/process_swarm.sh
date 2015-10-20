#!/bin/sh
#
# Run EEF chain on Swarm data

sats="A B C"

prog="./main"

prefix="log"

# set for profiles only
flags="-p"
#flags=""

for sat in ${sats}; do
  log_dir="${prefix}_${sat}"
  idx_file="swarm${sat}.idx"

  rm -rf ${log_dir}
  mkdir ${log_dir}
  mkdir ${log_dir}/plots

  echo "Processing satellite ${sat} (${idx_file}, ${log_dir})"
  screen -d -m -S swarm${sat} ${prog} -s ${idx_file} -l ${log_dir} ${flags}
done
