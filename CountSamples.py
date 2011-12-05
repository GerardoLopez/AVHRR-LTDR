#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import glob
import os
import sys

import numpy
import osgeo.gdal as gdal
from osgeo.gdalconst import *

def GetFileList(DataDir):
	FileList = glob.glob(DataDir + '/AVH09C1.A???????.N??.???.?????????????.tif')
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

#--------------------------------------------------------------------------------#
from IPython import embed

DataDir = sys.argv[1]

FileList, Year, DoY = GetFileList(DataDir)

# From the first file get dimensions
rows, cols, NumberOfBands = GetDimensions(FileList[0])

# Create aray where to store the NumberOfSamples
NumberOfSamples = numpy.zeros((rows, cols), numpy.float32)

FileNumber = 1
for File in FileList:
	print File
	dataset = gdal.Open(File, GA_ReadOnly)
	reflectance = dataset.GetRasterBand(1).ReadAsArray()

	indices = numpy.where(reflectance > 0)
	NumberOfSamples[indices] = NumberOfSamples[indices] + 1.0

	dataset = None

print "Writing results to a file..."
format = "GTiff"
driver = gdal.GetDriverByName(format)
new_dataset = driver.Create( 'NumberOfSamples.tif', cols, rows, GDT_Float32 )
new_dataset.GetRasterBand(1).WriteArray(NumberOfSamples)
new_dataset = None

