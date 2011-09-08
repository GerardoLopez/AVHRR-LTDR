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

	Year = numpy.zeros((len(FileList)), numpy.byte)
	DoY = numpy.zeros((len(FileList)), numpy.byte)

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

#--------------------------------------------------------------------------------#
DataDir = sys.argv[1]

FileList, Year, DoY = GetFileList(DataDir)

# From the first file get dimensions
rows, cols, NumberOfBands = GetDimensions(FileList[0])

# Create aray where to store the composite
composite = numpy.zeros((rows, cols), numpy.int16)

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
	#for band in range(1,2):
		FileNumber = 1
		for File in FileList:
			print File
			dataset = gdal.Open(File, GA_ReadOnly)
			print 0, InitRow, cols, (rows / NumberOfTiles)
			layerstack[:,:,band-1,FileNumber-1] = dataset.GetRasterBand(band).ReadAsArray(0, InitRow, cols, (rows / NumberOfTiles))

			FileNumber += 1

		format = "GTiff"
		driver = gdal.GetDriverByName(format)
		new_dataset = driver.Create( 'test_band' + str(band) +'.tif', cols, (rows / NumberOfTiles), NumberOfFiles, GDT_Int16 )

		print "Saving data band", band 
		for i in range(1,NumberOfFiles+1):
			new_dataset.GetRasterBand(i).WriteArray(layerstack[:,:,band-1,i-1])
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



