#!/bin/sh
#
# Make index file for CHAMP

datadir="$DATAHOME/CHAMP/Stage1_CHAOS"
idxfile="champ.idx"

find $datadir -name "*.cdf" | sort -g > ${idxfile}
