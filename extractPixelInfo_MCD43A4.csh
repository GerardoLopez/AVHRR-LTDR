#!/bin/csh

set FUNC = extractPixelInfo_MCD43A4.csh
set init_time = `date`

echo ""
echo ${FUNC}: Checking command line arguments

if ($#argv != 3) then
  echo ""
  echo Usage: $FUNC SAMPLE LINE OUTPUTFILE
  exit 1
endif

set sample = $1
set line = $2
set TextFile = $3

set SRCDIR = $HOME/LDRAG/src
set DATADIR = "/Users/glopez/LDRAG/data/MODIS/MCD43A4"

echo "Year,DoY,Nadir_Reflectance_Band1,Nadir_Reflectance_Band2,Nadir_Reflectance_Band3,Nadir_Reflectance_Band4,Nadir_Reflectance_Band5,Nadir_Reflectance_Band6,Nadir_Reflectance_Band7" >> $TextFile

foreach LTDR_FILE (`ls $DATADIR/????/*hdf`)
	echo $LTDR_FILE
	set filename = `basename $LTDR_FILE`
	set Year = `echo $filename | cut -d. -f2 | cut -c2-5`
	set DoY = `echo $filename | cut -d. -f2 | cut -c6-8`

	$SRCDIR/extractPixelInfo.py $LTDR_FILE $sample $line $sample $line $Year $DoY >> $TextFile
end
