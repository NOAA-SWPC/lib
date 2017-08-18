#!/bin/sh
#
# Run stage1 for a satellite and all years

dmsp_sat="F15"
if test -n "$1"; then
  dmsp_sat="$1"
fi

for year in `seq 2009 2015`; do
  screen -d -m -S p${dmsp_sat}_${year} sh stage1.sh ${dmsp_sat} ${year}
done
