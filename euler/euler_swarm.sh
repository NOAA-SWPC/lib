#!/bin/sh
#
# This script takes 3 separate Swarm-formatted euler angle files
# and combines them to a single file. Each of the 3 files must
# have exactly the same timestamps in the first column

outfile="MSW_EUL_BOUMME-0.1.txt"

# generate individual files
prog="./print"

fileA="MSWA_EUL_BOUMME-0.1.txt"
fileB="MSWB_EUL_BOUMME-0.1.txt"
fileC="MSWC_EUL_BOUMME-0.1.txt"

${prog} -i euler.0 -s $fileA
${prog} -i euler.1 -s $fileB
${prog} -i euler.2 -s $fileC

# store header first
cat $fileA | grep "^#" > $outfile

file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)

cat $fileA | grep -v "^#" | awk '{printf "%10.2f %11.6f %11.6f %11.6f\n", $1,$2,$3,$4}' > $file1
cat $fileB | grep -v "^#" | awk '{printf "%11.6f %11.6f %11.6f\n", $2,$3,$4}' > $file2
cat $fileC | grep -v "^#" | awk '{printf "%11.6f %11.6f %11.6f\n", $2,$3,$4}' > $file3

paste -d " " $file1 $file2 $file3 >> $outfile

rm -f $file1
rm -f $file2
rm -f $file3

echo "output file is $outfile"
