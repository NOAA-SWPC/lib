#!/bin/sh
#
# Run process_year on all years making a separate screen for each
# year

for year in $(seq 1970 1979); do
  screen -d -m -S p${year} sh process_year.sh ${year}
done
