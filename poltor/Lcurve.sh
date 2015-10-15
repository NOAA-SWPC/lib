#!/bin/sh
#
# L-curve analysis

data_file="data/champ_2000_2010_15s.dat"

# number of robust iterations
niter="1"

# starting/ending values for damping parameter
param_start="0.0"
param_end="0.1"

# number of steps to take
n="20"
param_step=$(echo "scale=4; (${param_end} - ${param_start}) / (${n} - 1)" | bc -l)

# output file
outfile="Lcurve.dat"

rm -f ${outfile}
echo "# Field 1: log ||b - A x||" >> ${outfile}
echo "# Field 2: log ||L x||" >> ${outfile}
echo "# Field 3: alpha_int" >> ${outfile}
echo "# Field 4: alpha_sh" >> ${outfile}
echo "# Field 5: alpha_tor" >> ${outfile}

for i in $(seq 1 ${n}); do
  alpha=$(echo "${param_start} + ${param_step}*($i-1)" | bc -l)
  echo "Running with alpha = ${alpha}"
  ./invert -i ${data_file} -q ${niter} --alpha_int ${alpha} --lcurve_file ${outfile} -l LLS.dat
done
