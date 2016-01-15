#!/bin/sh
# Create final Est/Ist data file from previously computed files

outfile="/data/Indices/ESTIST/Est_Ist_index.pli"

proc_script="/data/palken/lib/estist/process_year.sh"
datadir="/data/palken/lib/estist/data"

# Generate new data file for current year
sh ${proc_script} 2016

echo "Generating ${outfile}"
cat ${datadir}/Est_Ist*.dat | grep -v "\-9999" > ${outfile}
