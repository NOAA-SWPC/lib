#!/bin/sh
# Compute rms differences between MEA/DED and Weimer model for various storms

prog="./storm_rms"
tprog="../common/time2str"
outdir="output2"

rm -rf $outdir
mkdir $outdir

# October 30 2003
ts=1067472000
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_200310.dat -w ~/Downloads/mserv/MEA_Model_102003_IAGA2002.txt -t $ts > dat_${ts}_MEA 2> $outdir/output_200310_MEA

# November 20 2003
ts=1069286400
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_200310.dat -w ~/Downloads/mserv/MEA_Model_112003_IAGA2002.txt -t $ts > dat_${ts}_MEA 2> $outdir/output_200311_MEA

# March 17 2015 / MEA
ts=1426550400
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_032015_IAGA2002.txt -t $ts > dat_${ts}_MEA 2> $outdir/output_201503_MEA

# June 23 2015 / MEA
ts=1435017600
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_062015_IAGA2002.txt -t $ts > dat_${ts}_MEA 2> $outdir/output_201506_MEA

# Oct 7 2015 / MEA
ts=1444197600
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_102015_IAGA2002.txt -t $ts > dat_${ts}_MEA 2> $outdir/output_201510_MEA

# Dec 20 2015 / MEA
ts=1450594800
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_122015_IAGA2002.txt -t $ts > dat_${ts}_MEA 2> $outdir/output_20151220_MEA

# Dec 31 2015 / MEA
ts=1451545200
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/MEA/MEA_2015.dat -w ~/Downloads/mserv/MEA_Model_122015_IAGA2002.txt -t $ts > dat_${ts}_MEA 2> $outdir/output_20151231_MEA

###################################### DED ###########################

# March 17 2015 / DED
ts=1426550400
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_032015_IAGA2002.txt -t $ts > dat_${ts}_DED 2> $outdir/output_201503_DED

# June 23 2015 / DED
ts=1435017600
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_062015_IAGA2002.txt -t $ts > dat_${ts}_DED 2> $outdir/output_201506_DED

# Oct 7 2015 / DED
ts=1444197600
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_102015_IAGA2002.txt -t $ts > dat_${ts}_DED 2> $outdir/output_201510_DED

# Dec 20 2015 / DED
ts=1450594800
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_122015_IAGA2002.txt -t $ts > dat_${ts}_DED 2> $outdir/output_20151220_DED

# Dec 31 2015 / DED
ts=1451545200
echo "running date $($tprog $ts)"
$prog -i ~/Downloads/mserv/data/DED/DED_2015.dat -w ~/Downloads/mserv/DED_Model_122015_IAGA2002.txt -t $ts > dat_${ts}_DED 2> $outdir/output_20151231_DED

echo "Output dir is $outdir"
