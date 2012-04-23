#!/bin/csh

set FUNC = extractPixelInfo
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
set DATADIR = "/Users/glopez/LDRAG/data/MODIS/MCD15A2"

echo "Year,DoY,Fpar_1km,Lai_1km,FparLai_QC,FparExtra_QC,FparStdDev_1km,LaiStdDev_1km" >> $TextFile

foreach LTDR_FILE (`ls $DATADIR/????/*hdf`)

	echo $LTDR_FILE
	set filename = `basename $LTDR_FILE`
	set Year = `echo $filename | cut -d. -f2 | cut -c2-5`
	set DoY = `echo $filename | cut -d. -f2 | cut -c6-8`

	$SRCDIR/extractPixelInfo.py $LTDR_FILE $sample $line $sample $line $Year $DoY >> $TextFile
end
