#!/bin/sh
#
# Plot declination/inclination maps
#
# Usage: gmtplot_di.sh [xyzfile]

xyzfile="dat.xyz"
if test -n "$1"; then
  xyzfile="$1"
fi

cptfile="norm0.1.cpt"
#cptfile="norm0.05.cpt"
cbstep="0.05"
cbmax="0.1"

function dogmt
{
  proj="$1"
  region="$2"
  gflags="$3"

  # range of colorbar and step size in nT
  #cbrange="-10/10"

  # Convert xyz to grid file
  grdfile="${xyzfile}.grd"
  xyz2grd ${xyzfile} -G${grdfile} -R-180/180/-90/90 -I100+/100+

  # make color file
  cptfile="colors.cpt"
  makecpt -Cpolar -T-${cbmax}/${cbmax}/${cbstep} -Z > ${cptfile}

  # Make image
  grdimage -E300 ${gflags} -K -P -B0g0/0g0 ${grdfile} -C${cptfile} -J${proj} -R${region} >> $outfile
  pscoast -O -P -A10/1 -Di -W2 -B45g0f0/10g0f0 -R${region} -J${proj} -K >> $outfile
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
psscale -O -D-1i/-1i/5i/0.25ih -B${cbstep}/:deg: -C${cptfile} >> $outfile

echo "Output is ${outfile}"
