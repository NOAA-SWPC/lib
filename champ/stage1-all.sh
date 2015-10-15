#!/bin/sh
#
# Run stage1 on all data files making a separate screen for each
# year

for year in $(seq 2000 2010); do
  screen -d -m -S p${year} sh stage1.sh ${year}
done
