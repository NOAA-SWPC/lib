#!/bin/sh
#
# Prepare data files for plotting with plot_euler.gp

suffix="6am6pm_2"

function prepare_files
{
  time_period="$1"
  prefix="res_${time_period}_${suffix}"

  # Filter QD latitude to get LT/satdir of equator crossing
  outfile="res_${time_period}_qdfilt"
  echo "QD filtering ${prefix} files (output prefix: ${outfile})..."
  cat ${prefix}.sat0 | datasel -c 6 --min -5 --max 5 > ${outfile}.sat0
  cat ${prefix}.sat1 | datasel -c 6 --min -5 --max 5 > ${outfile}.sat1
  cat ${prefix}.sat2 | datasel -c 6 --min -5 --max 5 > ${outfile}.sat2
}

prepare_files "10d"
