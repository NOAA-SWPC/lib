#!/bin/sh
#
# Make index file for CHAMP

datadir="/nfs/satmag/CHAMP/Stage2"
idxfile="champ.idx"

find $datadir -name "*.cdf" | sort -g > ${idxfile}
