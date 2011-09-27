#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import glob
import os
import sys

import numpy
import osgeo.gdal as gdal
from osgeo.gdalconst import *

from scipy.linalg import lstsq, inv

def GetFileList(DataDir):
	FileList = glob.glob(DataDir + '/AVH09C1.A???????.N??.???.?????????????.img')
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

def BRDF_Inverter(Reflectance, Kernels, DoY, NumberOfParameters):
	rows, cols, NumberOfBands, days = Reflectance.shape
	parameters = numpy.zeros((rows, cols, NumberOfBands, NumberOfParameters), numpy.int16)
	#(900, 7200, 10, 31)
	for j in range(0,cols):
		for i in range(0,rows):
			# Extract temporal profile for SREFL CH1, CH2 and CH3
			ReflectanceProfile = Reflectance[i,j,:,:]
			# Based on valid observations on SREFL_CH1
			NumberOfSamples = len(numpy.where( ReflectanceProfile[0] > 0)[0])

			if NumberOfSamples >= 1:
				IndexData = numpy.where(ReflectanceProfile[0] > 0)[0]
				refl = ReflectanceProfile[:,IndexData]
				K = numpy.ones([3, NumberOfSamples])
				K[1,:] = Kernels[i,j,0,IndexData]
				K[2,:] = Kernels[i,j,1,IndexData]

				(P, rho_residuals, rank, svals) = lstsq ( K.T, refl.T)
				parameters[i,j,:,:] = P

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

# Create aray where to store the composite
# all bands in LTDR plus DoY of observation, NumberOfSamples
parameters = numpy.zeros((rows, cols, NumberOfBands, NumberOfParameters), numpy.int16)

# Depending on the processing system, the composite could be created storing ALL
# datasets in RAM, however for prototyping a tile-based processing will be implemented
# 3 tiles will be the default setting
NumberOfTiles = 3

for Tile in range(1,NumberOfTiles+1):
	InitRow = (Tile - 1) * (rows / NumberOfTiles)
	EndRow = (Tile * (rows / NumberOfTiles)) - 1
	print "Processing rows:", InitRow, "to", EndRow

	# Create tmp layerstack
	NumberOfFiles = len(FileList)
	Reflectance = numpy.zeros(((rows / NumberOfTiles), cols, NumberOfBands, NumberOfFiles), numpy.int16)
	Kernels = numpy.zeros(((rows / NumberOfTiles), cols, 2, NumberOfFiles), numpy.float32)

	FileNumber = 1
	for File in FileList:
		print File
		dataset = gdal.Open(File, GA_ReadOnly)
		for band in range(1, NumberOfBands+1):
			BandData = dataset.GetRasterBand(band).ReadAsArray()
			Reflectance[:,:,band-1,FileNumber-1] = BandData[InitRow:EndRow+1,0:cols]
			BandData = None

		dataset = None

		KernelsFile = File[0:len(File)-3] + "kernels.img"
		for band in range(1,3):
			dataset = gdal.Open(KernelsFile, GA_ReadOnly)
			BandData = dataset.GetRasterBand(band).ReadAsArray()
			Kernels[:,:,band-1,FileNumber-1] = BandData[InitRow:EndRow+1,0:cols]

		dataset = None

		FileNumber += 1

	print "Performing BRDF model inversion..."
	#ipshell = embed()

	tmpParameters = BRDF_Inverter(Reflectance, Kernels, DoY, NumberOfParameters)
	parameters[InitRow:EndRow+1,0:cols,:,:] = tmpParameters
	tmpParameters = None

print "Writing results to a file..."
format = "GTiff"
driver = gdal.GetDriverByName(format)
new_dataset = driver.Create( 'BRDF_Parameters.tif', cols, rows, NumberOfBands*NumberOfParameters, GDT_Int16 )

k = 1
for i in range(1,NumberOfBands):
	for j in range(1,NumberOfParameters):
		new_dataset.GetRasterBand(k).WriteArray(parameters[:,:,i-1,j-1])
		k += 1

new_dataset = None

