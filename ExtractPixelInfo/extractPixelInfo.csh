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

set SRCDIR = $HOME/GCII/src/AVHRR-LTDR/ExtractPixelInfo
set DATADIR = `ls -d /Users/glopez/GCII/data/LTDR/Filtered`

echo "DoY,SREFL_CH1,SREFL_CH2,SREFL_CH3,BT_CH3,BT_CH4,BT_CH5,SZEN,VZEN,RELAZ,QA" >> $TextFile

foreach LTDR_FILE (`ls -d $DATADIR/AVH09C1*.img`)

	echo $LTDR_FILE
	set filename = `basename $LTDR_FILE`
	set DoY = `echo $filename | cut -d. -f2 | cut -c6-8`

	$SRCDIR/extractPixelInfo.py $LTDR_FILE $sample $line $sample $line $DoY >> $TextFile

end
