#!/bin/bash

FUNC=ProcessExtractSDS_FromLTDR.sh
echo ""
echo ${FUNC}: Checking command line arguments

if [ $# -ne 3 ]; then
	echo ""
	echo Usage: $FUNC DATADIR INITDAY ENDDAY
	echo To process July, 1998 LTDR on /Users/glopez/GCII/data/LTDR/N14/1998
	echo $FUNC /Users/glopez/GCII/data/LTDR/N14/1998 183 213
	echo ""
	exit 1
fi

init_time=`date`
echo $init_time

SRCDIR=$HOME/GCII/src/AVHRR-LTDR

DATADIR=$1
INITDAY=$2
ENDDAY=$3

for SDS in `ls $DATADIR/*.hdf`; do
	DoY=`echo $SDS | cut -d. -f2 | cut -c6-8`
	if [ $DoY -ge $INITDAY ] && [ $DoY -le $ENDDAY ]; then
		echo "Processing $SDS..."
		echo $SRCDIR/ExtractSDS_FromLTDR.sh $SDS
		$SRCDIR/ExtractSDS_FromLTDR.sh $SDS
	fi

done

echo ""
echo " Started  at $init_time"
echo " Finished at "`date`
