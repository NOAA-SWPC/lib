#!/bin/sh
#
# Run process_year on all years making a separate screen for each
# year

for year in $(seq 1990 1999); do
  screen -d -m -S p${year} sh process_year.sh ${year}
done
