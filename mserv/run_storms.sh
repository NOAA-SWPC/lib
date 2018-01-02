#!/bin/sh
# Compute rms differences between MEA/DED and Weimer model for various storms

prog="./storm_rms"
tprog="time2str"
outdir="output"

rm -rf $outdir
mkdir $outdir

# October 30 2003
ts=1067472000
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_200310.dat -w ~/Downloads/mserv/MEA_Model_102003_IAGA2002.txt -t $ts -o $outdir/MEA_200310.txt 2> $outdir/MEA_200310.stderr

# November 20 2003
ts=1069286400
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_200310.dat -w ~/Downloads/mserv/MEA_Model_112003_IAGA2002.txt -t $ts -o $outdir/MEA_200311.txt 2> $outdir/MEA_200311.stderr

# March 17 2015 / MEA
ts=1426550400
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_032015_IAGA2002.txt -t $ts -o $outdir/MEA_201503.txt 2> $outdir/MEA_201503.stderr

# June 23 2015 / MEA
ts=1435017600
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_062015_IAGA2002.txt -t $ts -o $outdir/MEA_201506.txt 2> $outdir/MEA_201506.stderr

# Oct 7 2015 / MEA
ts=1444197600
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_102015_IAGA2002.txt -t $ts -o $outdir/MEA_201510.txt 2> $outdir/MEA_201510.stderr

# Dec 20 2015 / MEA
ts=1450594800
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_122015_IAGA2002.txt -t $ts -o $outdir/MEA_20151220.txt 2> $outdir/MEA_20151220.stderr

# Dec 31 2015 / MEA
ts=1451545200
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_122015_IAGA2002.txt -t $ts -o $outdir/MEA_20151231.txt 2> $outdir/MEA_20151231.stderr

###################################### DED ###########################

# March 17 2015 / DED
ts=1426550400
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_032015_IAGA2002.txt -t $ts -o $outdir/DED_201503.txt 2> $outdir/DED_201503.stderr

# June 23 2015 / DED
ts=1435017600
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_062015_IAGA2002.txt -t $ts -o $outdir/DED_201506.txt 2> $outdir/DED_201506.stderr

# Oct 7 2015 / DED
ts=1444197600
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_102015_IAGA2002.txt -t $ts -o $outdir/DED_201510.txt 2> $outdir/DED_201510.stderr

# Dec 20 2015 / DED
ts=1450594800
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_122015_IAGA2002.txt -t $ts -o $outdir/DED_20151220.txt 2> $outdir/DED_20151220.stderr

# Dec 31 2015 / DED
ts=1451545200
echo "running date $($tprog -t $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_122015_IAGA2002.txt -t $ts -o $outdir/DED_20151231.txt 2> $outdir/DED_20151231.stderr

echo "Output dir is $outdir"
