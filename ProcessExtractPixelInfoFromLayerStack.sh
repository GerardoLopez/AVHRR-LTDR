#!/bin/bash

DATADIR="/Users/glopez/GCII/data/ATSR-WFA"

for TimePeriod in `ls -d $DATADIR/*_*`
do
	echo $TimePeriod
	#cd $TimePeriod
	filename=`basename $TimePeriod`
	NumberOfFires=`wc -l $TimePeriod/${filename}.txt | awk '{print $1}'`
	for i in `seq 1 1 $NumberOfFires`
	do
		fire=`awk -v record=$i '{if (NR==record) print $0}' $TimePeriod/${filename}.txt`

		sample=`echo $fire | cut -d" " -f1`
		line=`echo $fire | cut -d" " -f2`
		NumberOfFires=`echo $fire | cut -d" " -f3`
		
		day=`echo $fire | cut -d" " -f4`
		month=`echo $fire | cut -d" " -f5`
		year=`echo $fire | cut -d" " -f6`

		JulianDay=`$HOME/GCII/src/AVHRR-LTDR/convert_gregorian_date_to_julian_day.csh $day $month $year`

		echo $fire
		$HOME/GCII/src/AVHRR-LTDR/extractPixelInfoFromLayerStack.py $sample $line $JulianDay $NumberOfFires
	done
done


