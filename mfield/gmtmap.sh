#!/bin/sh
#
# Plot X,Y,Z,F maps of Swarm residuals

prefix="res"

# Converted to EPS later
outfile="${prefix}.ps"

cptfile="$GMTHOME/norm10.cpt"
cbstep="5"

# Columns of residual data file
loncol="5"
latcol="6"
qdcol="7"
Fcol="8"
Xcol="9"
Ycol="10"
Zcol="11"
FMF="16"
XMF="17"
YMF="18"
ZMF="19"

# Args:
# projection
# region
# flags to pass to grdimage/pscoast
# plot title if any
# xyzfile
# zcol: component to plot
# title2
# mfcol: flag to check if data used for MF modeling
function dogmt
{
  proj="$1"
  region="$2"
  gflags="$3"
  title="$4"
  infile="$5"
  zcol="$6"
  title2="$7"
  mfcol="$8"

  if [ -n "${title}" ]; then
    pscoast_flags="-B:.${title}:"
  else
    pscoast_flags=""
  fi

  xyzfile=$(mktemp)
  cat ${infile} | datasel -c $mfcol --eq 1 | awk '{print $x,$y,$z}' x="${loncol}" y="${latcol}" z="${zcol}" > ${xyzfile}

  # Convert xyz to grid file
  grdfile="${xyzfile}.grd"
  gmt blockmean ${xyzfile} -R-180/180/-90/90 -I150m > ${xyzfile}.bm
  gmt xyz2grd ${xyzfile}.bm -G${grdfile} -R-180/180/-90/90 -I150m

  # Make image
  gmt grdimage -E300 ${gflags} -K -P -B0g0/0g0 ${grdfile} -C${cptfile} -J${proj} -R${region} >> $outfile
  #pscoast ${gflags} -O -P -A10/1 -Di -W2 -R${region} -J${proj} -K ${pscoast_flags} >> $outfile
  gmt pscoast ${gflags} -O -P -A5000 -Di -W2 -R${region} -J${proj} -K ${pscoast_flags} >> $outfile

  if [ -n "${title2}" ]; then
    gmt pstext ${gflags} -Xa-0.5 -Ya3.5 -N -R${region} -J${proj} -P -O -K << EOF >> $outfile
-180  100   14  0  1  RB  ${title2}
EOF
  fi

  rm -f ${grdfile}
  rm -f ${xyzfile}
  rm -f ${xyzfile}.bm
}

# select EPS file output
gmt set PS_MEDIA Custom_10ix8i

# Size of plot titles
gmt set FONT_TITLE 20

# Size of colorbar labels
gmt set FONT_ANNOT_PRIMARY 10

rm -f $outfile

# Equatorial view
p="W0/7"
r="-180/180/-90/90"

file="${prefix}.sat0"
echo "Rendering ${file}..."
dogmt $p $r "-Xa0.8i -Ya5.8i" "Swarm-A" ${file} ${Xcol} "" ${XMF}
dogmt $p $r "-Xa0.8i -Ya4.2i -O" "" ${file} ${Ycol} "" ${YMF}
dogmt $p $r "-Xa0.8i -Ya2.6i -O" "" ${file} ${Zcol} "" ${ZMF}
dogmt $p $r "-Xa0.8i -Ya1.0i -O" "" ${file} ${Fcol} "" ${FMF}

file="${prefix}.sat1"
echo "Rendering ${file}..."
dogmt $p $r "-Xa3.8i -Ya5.8i -O" "Swarm-B" ${file} ${Xcol} "" ${XMF}
dogmt $p $r "-Xa3.8i -Ya4.2i -O" "" ${file} ${Ycol} "" ${YMF}
dogmt $p $r "-Xa3.8i -Ya2.6i -O" "" ${file} ${Zcol} "" ${ZMF}
dogmt $p $r "-Xa3.8i -Ya1.0i -O" "" ${file} ${Fcol} "" ${FMF}

file="${prefix}.sat2"
echo "Rendering ${file}..."
dogmt $p $r "-Xa6.8i -Ya5.8i -O" "Swarm-C" ${file} ${Xcol} "" ${XMF}
dogmt $p $r "-Xa6.8i -Ya4.2i -O" "" ${file} ${Ycol} "" ${YMF}
dogmt $p $r "-Xa6.8i -Ya2.6i -O" "" ${file} ${Zcol} "" ${ZMF}
dogmt $p $r "-Xa6.8i -Ya1.0i -O" "" ${file} ${Fcol} "" ${FMF}

# Add row labels - the -R and -J options correspond to the PAPER_MEDIA
# size command above
gmt pstext -N -Ri0/7/0/8 -JX7i/8i -P -O -K << EOF >> $outfile
0.5  6.4   14  0  1  RB  X
0.5  4.8   14  0  1  RB  Y
0.5  3.2   14  0  1  RB  Z
0.5  1.6   14  0  1  RB  F
EOF

# Color scale bar
gmt psscale -O -D3.5i/0.8i/6i/0.20ih -B${cbstep}/:nT: -C${cptfile} >> $outfile

outfile_png=$(echo ${outfile} | sed -e 's/ps/png/')
outfile_eps=$(echo ${outfile} | sed -e 's/ps/eps/')

# Convert to EPS to get rid of whitespace around edges
#echo "Converting to EPS..."
ps2eps -f ${outfile}

# Convert to PNG
#echo "Converting to PNG..."
convert -density 130 -flatten ${outfile_eps} ${outfile_png}

# Delete PS/EPS files
rm -f ${outfile}
rm -f ${outfile_eps}

echo "Output is ${outfile_png}"
