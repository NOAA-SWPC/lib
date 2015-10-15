#!/bin/sh
#
# Usage: gmtplot.sh [xyzfile]

xyzfile="dat.xyz"
if test -n "$1"; then
  xyzfile="$1"
fi

cbmax="2"
cbstep="1"
cblabel="nT"
cptfile="norm${cbmax}.cpt"

# make constant scale color file
cptfile="colors.cpt"
makecpt -Cpolar -T-${cbmax}/${cbmax}/${cbstep} -Z > ${cptfile}

overlay=""

function dogmt
{
  xyz="$1"
  proj="$2"
  region="$3"
  pos="$4"
  title="$5"

  # Convert xyz to grid file
  grdfile="${xyz}.grd"
  bmfile="${xyz}.bm"

  blockmean ${xyz} -R-180/180/-90/90 -I30m > ${bmfile}
  #xyz2grd ${bmfile} -G${grdfile} -R-180/180/-90/90 -I240m/240m

  surface ${bmfile} -R-180/180/-90/90 -I300m -G${grdfile} -T0.25 -C0.1 -VL
  #sphinterpolate ${xyz} -G${grdfile} -R-180/180/-90/90 -I250m/250m -Q2 -T -V

  # Make image
  grdimage -E300 ${pos} ${overlay} -K -B0g0/0g0 ${grdfile} -C${cptfile} -J${proj} -R${region} >> $outfile
  pscoast ${pos} -O -A10/1 -Di -W2 -B45g0f0/10g0f0:."${title}": -R${region} -J${proj} -K >> $outfile

  if [ -z "${overlay}" ]; then
    overlay="-O"
  fi

  rm -f ${grdfile}
  rm -f ${bmfile}
}

function plotcomp
{
  xyz="$1"
  xpos="$2"
  lon="$3"
  lat="$4"
  comp="$5"
  title="$6"

  xyztmp=$(mktemp)
  cat ${xyz} | grep -v "^#" | awk -v lon=${lon} -v lat=${lat} -v comp=${comp} '{print $lon,$lat,$comp}' > ${xyztmp}

  # Equatorial view
  #p="W0/6.8"
  p="W0/3.8"
  r="-180/180/-90/90"
  dogmt ${xyztmp} $p $r "-Xa${xpos}i -Ya5.0i" "${title}"

  # North pole projection
  xpos_np=$(echo "${xpos} + 0.25" | bc)
  p="G0/90/1.3"
  r="-180/180/40/90"
  dogmt ${xyztmp} $p $r "-Xa${xpos_np}i -Ya3.0i" ""

  # South pole projection
  xpos_sp=$(echo "${xpos_np} + 2.0" | bc)
  p="G0/-90/1.3"
  r="-180/180/-90/-40"
  dogmt ${xyztmp} $p $r "-Xa${xpos_sp}i -Ya3.0i" ""

  rm -f ${xyztmp}
}

outfile="${xyzfile}.ps"
rm -f $outfile

gmtset ANNOT_FONT_SIZE 10p
gmtset PAPER_MEDIA legal

# X component
echo "plotting X..."
plotcomp ${xyzfile} 0.5 5 6 12 "X"

# Y component
#echo "plotting Y..."
#plotcomp ${xyzfile} 5.0 5 6 13 "Y"

# Z component
#echo "plotting Z..."
#plotcomp ${xyzfile} 9.5 5 6 14 "Z"

# Color scale bar
psscale -O -D7i/2i/7i/0.5ih -B${cbstep}/:${cblabel}: -C${cptfile} >> $outfile

echo "Output is ${outfile}"
