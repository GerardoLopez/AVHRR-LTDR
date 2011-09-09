#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import glob
import os
import sys

import numpy
import osgeo.gdal as gdal
from osgeo.gdalconst import *

def GetFileList(DataDir):
	FileList = glob.glob(DataDir + '/*.img')
	FileList.sort()

	Year = numpy.zeros((len(FileList)), numpy.int16)
	DoY = numpy.zeros((len(FileList)), numpy.int16)

	i = 0
	for File in FileList:
		# Get Year and DoY from filename
		YearOfObservation = os.path.basename(File).split('.')[1][1:5]
		DoYOfObservation = os.path.basename(File).split('.')[1][5:8]

		Year[i] = YearOfObservation
		DoY[i] = DoYOfObservation

		i += 1

	return FileList, Year, DoY


def GetDimensions(File):
	dataset = gdal.Open(File, GA_ReadOnly)
	# Usually the AVH09C1 layerstack should contain 10 bands
	rows, cols, NumberOfBands = dataset.RasterYSize, dataset.RasterXSize, dataset.RasterCount
	dataset = None

	return rows, cols, NumberOfBands

def CreateComposite(layerstack, DoY):
	rows, cols, NumberOfBands, days = layerstack.shape
	composite = numpy.zeros((rows, cols, NumberOfBands+2), numpy.int16)
	#(900, 7200, 10, 31)
	for j in range(0,cols):
		for i in range(0,rows):
			# Extract temporal profile for SREFL_CH1 and SREFL_CH2, bands 1, 2 and 5
			profile_SREFL_CH1 = layerstack[i,j,0,:]
			profile_SREFL_CH2 = layerstack[i,j,1,:]
			profile_BT_CH4 = layerstack[i,j,4,:]

			NumberOfSamples = len(numpy.where(profile_SREFL_CH1 > 0)[0])

			if NumberOfSamples == 1:
				IndexData = numpy.where(profile_SREFL_CH1 > 0)[0]
				for band in range(0,NumberOfBands):
					composite[i,j,band] =  layerstack[i,j,band,IndexData]

				composite[i,j,NumberOfBands] = DoY[IndexData]
				composite[i,j,NumberOfBands+1] = NumberOfSamples

			IndexData = numpy.where(profile_SREFL_CH1 > 0)[0]
			NewProfile_SREFL_CH1 = profile_SREFL_CH1[IndexData]
			NewProfile_SREFL_CH2 = profile_SREFL_CH2[IndexData]
			NewProfile_BT_CH4 = profile_BT_CH4[IndexData]
			NewDoY = DoY[IndexData]
			# Get SREFL_CH1 + SREFL_CH2
			Sum_SREFL_CH_1_2 = NewProfile_SREFL_CH1 + NewProfile_SREFL_CH2

			if NumberOfSamples == 2:
				# Get the index of the lowest sum
				IndexLowestSum = numpy.argmin(Sum_SREFL_CH_1_2)

				for band in range(0,NumberOfBands):
					layerstack[i,j,band,IndexData[IndexLowestSum]]

				composite[i,j,NumberOfBands] = DoY[IndexData[IndexLowestSum]]
				composite[i,j,NumberOfBands+1] = NumberOfSamples

			if NumberOfSamples >= 3:
				# Get the indices of the 3 lowest sums
				NumberOfSamplesLowestSums = 3
				IndicesLowestSum = numpy.argsort(Sum_SREFL_CH_1_2)[0:NumberOfSamplesLowestSums]
				# From the above samples, get the index of the one with lowest BT_CH4
				IndicesLowestSumSorted = numpy.sort(IndicesLowestSum)
				IndexLowest_BT_CH4 = numpy.argmin(NewProfile_BT_CH4[IndicesLowestSumSorted])

				for band in range(0,NumberOfBands):
					layerstack[i,j,band,IndexData[IndicesLowestSumSorted][IndexLowest_BT_CH4]]

				composite[i,j,NumberOfBands] = DoY[IndexData[IndicesLowestSumSorted][IndexLowest_BT_CH4]]
				composite[i,j,NumberOfBands+1] = NumberOfSamples

				#ipshell = embed()

	return composite


#--------------------------------------------------------------------------------#
from IPython import embed

DataDir = sys.argv[1]

FileList, Year, DoY = GetFileList(DataDir)

# From the first file get dimensions
rows, cols, NumberOfBands = GetDimensions(FileList[0])

# Create aray where to store the composite
# all bands in LTDR plus DoY of observation, NumberOfSamples
composite = numpy.zeros((rows, cols, NumberOfBands+2), numpy.int16)

# Depending on the processing system, the composite could be created storing ALL
# datasets in RAM, however for prototyping a tile-based processing will be implemented
# 4 tiles will be the default setting
NumberOfTiles = 4

for Tile in range(1,NumberOfTiles+1):
	InitRow = (Tile - 1) * (rows / NumberOfTiles)
	EndRow = (Tile * (rows / NumberOfTiles)) - 1
	print "Processing rows:", InitRow, "to", EndRow

	# Create tmp layerstack
	NumberOfFiles = len(FileList)
	layerstack = numpy.zeros(((rows / NumberOfTiles), cols, NumberOfBands, NumberOfFiles), numpy.int16)

	for band in range(1, NumberOfBands+1):
		print "Extracting band", band
	#for band in range(1,2):
		FileNumber = 1
		for File in FileList:
			#print File
			dataset = gdal.Open(File, GA_ReadOnly)
			#print 0, InitRow, cols, (rows / NumberOfTiles)
			layerstack[:,:,band-1,FileNumber-1] = dataset.GetRasterBand(band).ReadAsArray(0, InitRow, cols, (rows / NumberOfTiles))

			FileNumber += 1
	
	test = CreateComposite(layerstack, DoY)

	format = "GTiff"
	driver = gdal.GetDriverByName(format)
	new_dataset = driver.Create( 'composite.tif', cols, (rows / NumberOfTiles), NumberOfBands, GDT_Int16 )

	for i in range(1,NumberOfBands+1):
		new_dataset.GetRasterBand(i).WriteArray(composite[:,:,i-1])
	new_dataset = None

	exit()


layerstack = numpy.zeros((rows, cols, NumberOfBands), numpy.int16)

print "Reading bands..."
for i in range(1,NumberOfBands+1):
	layerstack[:,:,i-1] = dataset.GetRasterBand(i).ReadAsArray()

dataset = None

# QA bit-flag band is stored in band number 10 
QA = layerstack[:,:,9]
QA_binary = numpy.zeros((rows, cols), numpy.int16)
QA_flags = [128, 16512]

print "Screening data using QA flags..."
for flag in QA_flags:
	QA_binary = numpy.where(QA == flag, 1, QA_binary)


for i in range(1,NumberOfBands+1):
	layerstack[:,:,i-1] = layerstack[:,:,i-1] * QA_binary


# Save results to a ENVI binary file
format = "ENVI"
driver = gdal.GetDriverByName(format)
new_dataset = driver.Create( 'SDS_layerstack_masked.img', cols, rows, NumberOfBands, GDT_Int16 )

print "Saving data..."
for i in range(1,NumberOfBands+1):
	new_dataset.GetRasterBand(i).WriteArray(layerstack[:,:,i-1])

new_dataset = None



