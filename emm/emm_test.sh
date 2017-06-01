#!/bin/sh
#
# Generate EMM test output files

prog="./emm_test"
outdir="./emm_test_data"

# EMM coefficients in msynth format
coefdir="./coef2015"

start_year="2000"
end_year="2015"

rm -rf ${outdir}
mkdir -p ${outdir}

for year in $(seq ${start_year} ${end_year}); do
  outfile="${outdir}/emm_test_${year}.out"
  echo "generating ${outfile}..."
  ${prog} -e ${year} ${coefdir}/*.txt > ${outfile} 2>/dev/null
done
