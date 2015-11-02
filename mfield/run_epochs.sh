#!/bin/sh
#
# Run magnetic field model on all epochs

start_epoch="2010.0"
end_epoch="2014.0"
numit="5"

prog="./mfield"
coef_dir="/nfs/satmag_work/palken/corr/mfield/coef_new"

lambda_sv="0.006"
lambda_sa="0.006"

# Data interval available for DMSP
data_start="2009.0"
data_end="2014.33"

rm -rf ${coef_dir}
mkdir ${coef_dir}

# step size = 1 month = 1/12 years
step="0.083333"

for epoch in $(seq ${start_epoch} ${step} ${end_epoch}); do

  diff1=$(echo ${epoch} - ${data_start} | bc)
  diff2=$(echo ${data_end} - ${epoch} | bc)

  # compute min(epoch-data_start, data_end-epoch, 1.5)

  # min(epoch-data_start,1.5)
  if [ 1 -eq $(echo "${diff1} < 1.5" | bc) ]; then
    min1=${diff1}
  else
    min1=1.5
  fi

  # min(min1, data_end-epoch)
  if [ 1 -eq $(echo "${min1} < ${diff2}" | bc) ]; then
    window_size=${min1}
  else
    window_size=${diff2}
  fi

  # find 3-year window centered on epoch
  low=$(echo ${epoch} - ${window_size} | bc)
  high=$(echo ${epoch} + ${window_size} | bc)

  echo "epoch = ${epoch} low = ${low} high = ${high} window_size = ${window_size}"

  # find start year and month
  start_year=${low/\.*}
  start_month=$(echo "(${low} - ${start_year}) * 12" | bc)
  start_month=$(printf "%.0f" ${start_month})

  end_year=${high/\.*}
  end_month=$(echo "(${high} - ${end_year}) * 12" | bc)
  end_month=$(printf "%.0f" ${end_month})

  # create index files of DMSP data
  echo "building index files for epoch ${epoch} (${start_month}/${start_year} to ${end_month}/${end_year})..."
  sh domodel.sh ${start_month} ${start_year} ${end_month} ${end_year} > /dev/null 2>&1

  coef_file="${coef_dir}/coef.${epoch}.dat"
  echo "computing model for epoch ${epoch}; output in ${coef_file}..."
  ${prog} -n ${numit} -e ${epoch} -c ${coef_file} -v ${lambda_sv} -a ${lambda_sa}
done
