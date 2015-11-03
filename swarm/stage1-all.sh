#!/bin/sh
#
# Run stage1 on all data files for one satellite, making a
# separate screen for each year
#
# Usage: stage1-all.sh [A|B|C]

sat="A"
if test -n "$1"; then
  sat="$1"
fi

for year in $(seq 2013 2015); do
  screen -d -m -S pSwarm${sat}_${year} sh stage1_chaos.sh ${sat} ${year}
done
