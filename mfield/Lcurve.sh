#!/bin/sh
#
# Attempt to find optimal damping parameter for SV/SA coeffs
# by plotting log ||f|| vs log ||x|| for different values of
# lambda

# model epoch
epoch="2014.4"

# SV damping parameter (fixed)
lambda_sv="0.0"

# number of robust iterations
niter="5"

# starting/ending values for lambda_sa
start="0.005"
end="0.01"

# number of steps to take
n="20"

step=$(echo "scale=4; (${end} - ${start}) / (${n} - 1)" | bc -l)

# output file
outfile="lambda.dat"

rm -f ${outfile}

echo "# Field 1: log ||f||" >> ${outfile}
echo "# Field 2: log ||x||" >> ${outfile}
echo "# Field 3: lambda_sv" >> ${outfile}
echo "# Field 3: lambda_sa" >> ${outfile}

for i in $(seq 1 ${n}); do
  #lambda_sa=$(echo "${start} + ${step}*($i-1)" | bc -l)
  #echo "Running with lambda_sa = ${lambda_sa}"
  #./mfield -e ${epoch} -n ${niter} -v ${lambda_sv} -a ${lambda_sa} >> ${outfile}

  lambda_sv=$(echo "${start} + ${step}*($i-1)" | bc -l)
  echo "Running with lambda_sv = ${lambda_sv}"
  ./mfield -e ${epoch} -n ${niter} -v ${lambda_sv} -l ${outfile}

  # Save spectrum for this lambda
  cp mfield.s mfield.s.${lambda_sv}
done
