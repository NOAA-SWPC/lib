#!/bin/sh
#
# Usage: gmtplot.sh [xyzfile]

xyzfile="dat.xyz"
if test -n "$1"; then
  xyzfile="$1"
fi

cbmax="20"
cbstep="5"
cptfile="norm${cbmax}.cpt"

# make constant scale color file
#cptfile="colors.cpt"
#makecpt -Cpolar -T-${cbmax}/${cbmax}/${cbstep} -Z > ${cptfile}

function dogmt
{
  proj="$1"
  region="$2"
  gflags="$3"

  xyztmp=$(mktemp)
  cat ${xyzfile} | grep -v "^#" | awk '{print $1,$2,$4}' > ${xyztmp}

  # Convert xyz to grid file
  grdfile="${xyzfile}.grd"

  blockmean ${xyztmp} -R-180/180/-90/90 -I30m > ${xyzfile}.bm
  xyz2grd ${xyzfile}.bm -G${grdfile} -R-180/180/-90/90 -I100+/100+

  #surface ${xyzfile}.bm -R-180/180/-90/90 -I300m -G${grdfile} -T0.25 -C0.1 -VL
  #sphinterpolate ${xyztmp} -G${grdfile} -R-180/180/-90/90 -I250m/250m -Q2 -T -V

  # Make image
  grdimage -E300 ${gflags} -K -P -B0g0/0g0 ${grdfile} -C${cptfile} -J${proj} -R${region} >> $outfile
  pscoast -O -P -A10/1 -Di -W2 -B45g0f0/10g0f0 -R${region} -J${proj} -K >> $outfile

  rm -f ${xyztmp}
  rm -f ${grdfile}
  rm -f ${xyzfile}.bm
}

outfile="${xyzfile}.ps"
rm -f $outfile

# Equatorial view
p="W0/6.8"
r="-180/180/-90/90"
dogmt $p $r "-X0.8i -Y6.5i"

# North pole projection
p="G0/90/2"
r="-180/180/40/90"
dogmt $p $r "-X0.4i -Y-2.6i -O"

# South pole projection
p="G0/-90/2"
r="-180/180/-90/-40"
dogmt $p $r "-X4.0i -Y0 -O"

# Color scale bar
psscale -O -D-1i/-1i/5i/0.25ih -B${cbstep}/:nT: -C${cptfile} >> $outfile

echo "Output is ${outfile}"
