#!/bin/sh
#
# This script converts the ASCII DMSP files to CDF format, additionally
# correcting the ephemeris from Bowman or TENA
#
# Usage: ./stage0.sh [dmsp_sat] [year]
#
# Errors are recorded to stage0.log file

# satellite to process
dmsp_sat="F15"

# year to process
year="2009"

if test -n "$1"; then
  dmsp_sat="$1"
fi

if test -n "$2"; then
  year="$2"
fi

echo "==== PROCESSING DMSP ${dmsp_sat} YEAR ${year} ===="

logdir="$MYLIBHOME/dmsp/log"
logfile="$logdir/stage0.${dmsp_sat}.${year}.log"

ascdir="$DATAHOME/DMSP/MAG/ASCII"
cdfdir="$DATAHOME/DMSP/MAG/Stage0"
prog="$MYLIBHOME/dmsp/stage0"

mkdir -p $cdfdir/$year
mkdir -p $logdir
rm -f ${logfile}

bowman()
{
  infile="$1"
  outfile="$2"

  # Bowman ephemeris data files
  bow_dir="/data/palken/DMSP/EPH_bowman"

  # find doy from inside DMSP file
  doy=`cat ${file} | awk 'NR == 7 {print $1}' | cut -c 5-7`

  bow_eph_file="${bow_dir}/${dmsp_sat}_${year}_ephemeris.gz"
  if [ ! -e "${bow_eph_file}" ]; then
    echo "===ERROR===: Bowman ephemeris file ${bow_eph_file} does not exist" >> ${logfile}
    return 1
  fi

  ${prog} -i ${file} -o ${outfile} -b ${bow_eph_file}
}

tena()
{
  infile="$1"
  outfile="$2"

  # Bowman ephemeris data files
  tena_dir="$DATAHOME/DMSP/EPH"

  # find doy from inside DMSP file
  doy=$(zcat ${infile} | awk 'NR == 7 {print $1}' | cut -c 5-7)

  tena_eph_file="${tena_dir}/${dmsp_sat}/${year}DMSP/D.${doy}.${year}.DOP"
  if [ ! -e "${tena_eph_file}" ]; then
    echo "===ERROR===: TENA ephemeris file ${tena_eph_file} does not exist" >> ${logfile}
    return 1
  fi

  ${prog} -i ${file} -o ${outfile} -t ${tena_eph_file}
}

# No ephemeris correction
orig()
{
  infile="$1"
  outfile="$2"

  ${prog} -i ${file} -o ${outfile}
}

cd $ascdir/$year
for file in `ls *${dmsp_sat}*`; do
  datestr=$(echo ${file} | sed -e "s/.*${year}\([0-9]\)\([0-9]\)\([0-9]\)\([0-9]\).*/${year}\1\2\3\4/g")
  outfile="${cdfdir}/${year}/DMSP_${dmsp_sat}_${datestr}_MFR_Stage0.cdf"

  # Check first if we already processed this file
  if [[ ! -f "${outfile}" ]]; then
    echo "Processing: ${file}"

    # Bowman ephemeris
    bowman ${file} ${outfile}

    # TENA ephemeris
    #tena ${file} ${outfile}

    # Original ephemeris
    #orig ${file} ${outfile}
  fi
done
