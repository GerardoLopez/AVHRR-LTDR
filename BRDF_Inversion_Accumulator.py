#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import glob
import os
import sys

import numpy
import osgeo.gdal as gdal
from osgeo.gdalconst import *

from scipy.linalg import lstsq, inv

def GetFileList(DataDir):
	FileList = glob.glob(DataDir + '/AVH09C1.A?????[2-3]?.N??.???.?????????????.tif')
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

def BRDF_Inverter(M, V, DoY, NumberOfParameters, NumberOfSamples):
	rows, cols, NumberOfBands = M.shape
	parameters = numpy.zeros((rows, cols, NumberOfBands, NumberOfParameters), numpy.float32)

	#(900, 7200, 10, 31)
	for j in range(0,cols):
		if j % 100 == 0:
			print "Processing columns", j
		for i in range(0,rows):
			if NumberOfSamples[i,j] >= 2:
				(P, rho_residuals, rank, svals) = lstsq(M, V)
				parameters[i,j,:,:] = P

				#if i == 925 and j == 1646:
				#ipshell = embed()

	return parameters


#--------------------------------------------------------------------------------#
from IPython import embed

DataDir = sys.argv[1]

FileList, Year, DoY = GetFileList(DataDir)

# From the first file get dimensions
rows, cols, NumberOfBands = GetDimensions(FileList[0])
# Overwrite number of bands with 3, only:
# SREFL_CH1, SREFL_CH2, SREFL_CH3
NumberOfBands = 3
NumberOfParameters = 3

ScaleFactor = 0.0001

# Create aray where to store the composite
# all bands in LTDR plus DoY of observation, NumberOfSamples
parameters = numpy.zeros((rows, cols, NumberOfBands, NumberOfParameters), numpy.float32)
NumberOfSamples = numpy.zeros((rows, cols), numpy.int16)

# Depending on the processing system, the composite could be created storing ALL
# datasets in RAM, however for prototyping a tile-based processing will be implemented
# 4 tiles will be the default setting
NumberOfTiles = 1

for Tile in range(1,NumberOfTiles+1):
	InitRow = (Tile - 1) * (rows / NumberOfTiles)
	EndRow = (Tile * (rows / NumberOfTiles)) - 1
	print "Processing rows:", InitRow, "to", EndRow

	# Create tmp layerstack
	NumberOfFiles = len(FileList)
	Reflectance = numpy.zeros(((rows / NumberOfTiles), cols, NumberOfBands), numpy.float32)
	Kernels = numpy.zeros(((rows / NumberOfTiles), cols, 2), numpy.float32)
	tmpNumberOfSamples = numpy.zeros((rows / NumberOfTiles, cols), numpy.int16)

	M = numpy.zeros(((rows / NumberOfTiles), cols, 3), numpy.float32)
	V = numpy.zeros(((rows / NumberOfTiles), cols, 3), numpy.float32)

	for File in FileList:
		print File
		dataset = gdal.Open(File, GA_ReadOnly)
		for band in range(1, NumberOfBands+1):
			BandData = dataset.GetRasterBand(band).ReadAsArray()
			if NumberOfTiles > 1:
				Reflectance[:,:,band-1] = BandData[InitRow:EndRow+1,0:cols] * ScaleFactor
			else:
				Reflectance[:,:,band-1] = BandData[:,:] * ScaleFactor

			# Number of samples is based on number of observations on band 2
			if band == 2:
				tmpNumberOfSamples[:,:] = numpy.where(BandData[:,:] > 0.0, tmpNumberOfSamples[:,:] + 1, tmpNumberOfSamples[:,:])

			BandData = None
		dataset = None

		KernelsFile = File[0:len(File)-3] + "kernels.tif"
		for band in range(1,3):
			dataset = gdal.Open(KernelsFile, GA_ReadOnly)
			BandData = dataset.GetRasterBand(band).ReadAsArray()
			if NumberOfTiles > 1:
				Kernels[:,:,band-1] = BandData[InitRow:EndRow+1,0:cols]
			else:
				Kernels[:,:,band-1] = BandData[:,:]

		dataset = None

		for j in range(0,cols):
			for i in range(0,rows):
				if tmpNumberOfSamples[i,j] >= 2:
					K = numpy.ones(3)
					K[1] = Kernels[i,j,0]
					K[2] = Kernels[i,j,1]

					R = Reflectance[i,j,:]

					M[i,j,:] = M[i,j,:] + (K.T * K)
					V[i,j,:] = V[i,j,:] + (K.T * R)

	print "Performing BRDF model inversion..."
	#ipshell = embed()

	#tmpParameters = BRDF_Inverter(Reflectance, Kernels, DoY, NumberOfParameters, tmpNumberOfSamples)
	tmpParameters = BRDF_Inverter(M, V, DoY, NumberOfParameters, tmpNumberOfSamples)
	parameters[InitRow:EndRow+1,0:cols,:,:] = tmpParameters
	NumberOfSamples[InitRow:EndRow+1,0:cols] = tmpNumberOfSamples
	tmpParameters = None
	tmpNumberOfSamples = None

print "Writing results to a file..."
format = "GTiff"
driver = gdal.GetDriverByName(format)
new_dataset = driver.Create( 'BRDF_Parameters_Accumulator.tif', cols, rows, (NumberOfBands*NumberOfParameters) + 1, GDT_Float32 )

k = 1
for i in range(1,NumberOfBands+1):
	for j in range(1,NumberOfParameters+1):
		new_dataset.GetRasterBand(k).WriteArray(parameters[:,:,j-1,i-1])
		k += 1

new_dataset.GetRasterBand((NumberOfBands*NumberOfParameters) + 1).WriteArray(NumberOfSamples[:,:])
new_dataset = None

