#!/bin/sh

coefdir="./coef"

prog="./mfield_dBdt"
outfile="dBdt_kou.dat"

rm -f ${outfile}

echo "# Field 1: time (years)" >> ${outfile}
echo "# Field 2: dX/dt (nT/year)" >> ${outfile}
echo "# Field 3: dY/dt (nT/year)" >> ${outfile}
echo "# Field 4: dZ/dt (nT/year)" >> ${outfile}

for file in $(ls ${coefdir}/*.dat); do
  ${prog} -c ${file} >> ${outfile}
done

echo "output file is ${outfile}"
