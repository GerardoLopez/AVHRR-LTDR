#!/bin/bash

SRCDIR=$HOME/GCII/src/AVHRR-LTDR
METADATADIR=$HOME/GCII/data/LTDR/metadata

LTDR_FILE=$1

# Extract SDS names from LTDR file
SDS_VARS=(`hdp dumpsds -h $1 | grep "Variable Name" | awk -F" = " '{print $2}'`)
NumberOfSDS=${#SDS_VARS[@]}
let NumberOfSDS=$NumberOfSDS-1

for i in `jot - 0 $NumberOfSDS`
do
	SDS=${SDS_VARS[$i]}

	# Export SDS to a generic binary file
	echo "hdp dumpsds -n $SDS -o $SDS.img -b $LTDR_FILE"
	hdp dumpsds -n $SDS -o $SDS.img -b $LTDR_FILE
	# Create ENVI header from template
	cp $METADATADIR/GenericHeader.hdr $SDS.hdr
done

# Create layersatck with all bands
BandList="SREFL_CH1.img SREFL_CH2.img SREFL_CH3.img BT_CH3.img BT_CH4.img BT_CH5.img SZEN.img VZEN.img RELAZ.img QA.img"
HeaderList="SREFL_CH1.hdr SREFL_CH2.hdr SREFL_CH3.hdr BT_CH3.hdr BT_CH4.hdr BT_CH5.hdr SZEN.hdr VZEN.hdr RELAZ.hdr QA.hdr"
echo "Creating temporal layerstack..."
gdal_merge.py -o SDS_layerstack.img -of ENVI -separate $BandList
rm $BandList $HeaderList

# Filter dataset using QA values
$SRCDIR/Filter_AVH09C1_by_QA.py
mv SDS_layerstack_masked.tif `basename $LTDR_FILE hdf`tif
#cp $METADATADIR/GenericHeaderLayerStack.hdr `basename $LTDR_FILE hdf`hdr

mv SDS_layerstack_kernels_masked.tif `basename $LTDR_FILE hdf`kernels.tif
#cp $METADATADIR/GenericHeaderKernels.hdr `basename $LTDR_FILE hdf`kernels.hdr

rm SDS_layerstack.??? SDS_layerstack_masked.hdr SDS_layerstack_kernels_masked.hdr

