#!/bin/sh
#
# Plot X,Y,Z,F maps of Swarm residuals

prefix="res_emm.dat"

# Converted to EPS later
outfile="${prefix}.ps"

cptfile="norm10.cpt"
cbstep="5"

# Args:
# projection
# region
# flags to pass to grdimage/pscoast
# plot title if any
# xyzfile
function dogmt
{
  proj="$1"
  region="$2"
  gflags="$3"
  title="$4"
  infile="$5"
  zcol="$6"
  title2="$7"
  vec="$8"

  if [ -n "${title}" ]; then
    pscoast_flags="-B:.${title}:"
  else
    pscoast_flags=""
  fi

  xyzfile=$(mktemp)
  if [ ${vec} -gt "0" ]; then
    cat ${infile} | datasel -c 6 --min -55 --max 55 | awk '{print $4,$5,$x}' x="${zcol}" > ${xyzfile}
  else
    # For some reason, there are large residuals of > 500 nT leaking
    # into dataset
    cat ${infile} | datasel -c ${zcol} --min -100 --max 100 | awk '{print $4,$5,$x}' x="${zcol}" > ${xyzfile}
  fi

  # Convert xyz to grid file
  grdfile="${xyzfile}.grd"
  blockmean ${xyzfile} -R-180/180/-90/90 -I150m > ${xyzfile}.bm
  xyz2grd ${xyzfile}.bm -G${grdfile} -R-180/180/-90/90 -I150m

  # Make image
  grdimage -E300 ${gflags} -K -P -B0g0/0g0 ${grdfile} -C${cptfile} -J${proj} -R${region} >> $outfile
  #pscoast ${gflags} -O -P -A10/1 -Di -W2 -R${region} -J${proj} -K ${pscoast_flags} >> $outfile
  pscoast ${gflags} -O -P -A5000 -Di -W2 -R${region} -J${proj} -K ${pscoast_flags} >> $outfile

  if [ -n "${title2}" ]; then
    pstext ${gflags} -Xa-0.5 -Ya3.5 -N -R${region} -J${proj} -P -O -K << EOF >> $outfile
-180  100   14  0  1  RB  ${title2}
EOF
  fi

  rm -f ${grdfile}
  rm -f ${xyzfile}
  rm -f ${xyzfile}.bm
}

# select EPS file output
gmtset PAPER_MEDIA Custom_7ix8i

# Size of plot titles
gmtset HEADER_FONT_SIZE 20

# Size of colorbar labels
gmtset ANNOT_FONT_SIZE_PRIMARY 10

rm -f $outfile

# Equatorial view
p="W0/1.8"
r="-180/180/-90/90"

file="${prefix}.sat0"
echo "Rendering ${file}..."
dogmt $p $r "-Xa0.8i -Ya4.6i" "Swarm-A" ${file} 8 "" 1
dogmt $p $r "-Xa0.8i -Ya3.4i -O" "" ${file} 9 "" 1
dogmt $p $r "-Xa0.8i -Ya2.2i -O" "" ${file} 10 "" 1
dogmt $p $r "-Xa0.8i -Ya1.0i -O" "" ${file} 7 "" 0

file="${prefix}.sat1"
echo "Rendering ${file}..."
dogmt $p $r "-Xa2.8i -Ya4.6i -O" "Swarm-B" ${file} 8 "" 1
dogmt $p $r "-Xa2.8i -Ya3.4i -O" "" ${file} 9 "" 1
dogmt $p $r "-Xa2.8i -Ya2.2i -O" "" ${file} 10 "" 1
dogmt $p $r "-Xa2.8i -Ya1.0i -O" "" ${file} 7 "" 0

file="${prefix}.sat2"
echo "Rendering ${file}..."
dogmt $p $r "-Xa4.8i -Ya4.6i -O" "Swarm-C" ${file} 8 "" 1
dogmt $p $r "-Xa4.8i -Ya3.4i -O" "" ${file} 9 "" 1
dogmt $p $r "-Xa4.8i -Ya2.2i -O" "" ${file} 10 "" 1
dogmt $p $r "-Xa4.8i -Ya1.0i -O" "" ${file} 7 "" 0

# Add row labels - the -R and -J options correspond to the PAPER_MEDIA
# size command above
pstext -N -Ri0/7/0/8 -JX7i/8i -P -O -K << EOF >> $outfile
0.5  5.0   14  0  1  RB  X
0.5  3.8   14  0  1  RB  Y
0.5  2.6   14  0  1  RB  Z
0.5  1.4   14  0  1  RB  F
EOF

# Color scale bar
psscale -O -D3.5i/0.8i/6i/0.20ih -B${cbstep}/:nT: -C${cptfile} >> $outfile

outfile_png=$(echo ${outfile} | sed -e 's/ps/png/')
outfile_eps=$(echo ${outfile} | sed -e 's/ps/eps/')

# Convert to EPS to get rid of whitespace around edges
echo "Converting to EPS..."
ps2eps -f ${outfile}

# Convert to PNG
echo "Converting to PNG..."
convert -density 130 -flatten ${outfile_eps} ${outfile_png}

# Delete PS/EPS files
rm -f ${outfile}
rm -f ${outfile_eps}

echo "Output is ${outfile_png}"
