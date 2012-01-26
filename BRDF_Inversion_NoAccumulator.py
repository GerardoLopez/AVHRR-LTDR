#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import glob
import os
import sys

import numpy
import math
import osgeo.gdal as gdal
from osgeo.gdalconst import *

from scipy.linalg import lstsq, inv

def NormalDistribution(mu,sigma):
	def f(x):
		z = 1.0*(x-mu)/sigma
		e = math.e**(-0.5*z**2)
		C = math.sqrt(2*math.pi)*sigma
		return 1.0*e/C
	return f


def GetPrior(PriorDataDir, strDoY, InitRow, EndRow):
	PriorFile = PriorDataDir + '/MCD43C2.Prior.' + strDoY + '.img'
	rows, cols, NumberOfBands = GetDimensions(PriorFile)

	rows = (EndRow-InitRow)+1
	Prior = numpy.zeros((rows,cols,NumberOfBands), numpy.float32)

	dataset = gdal.Open(PriorFile, GA_ReadOnly)

	for band in range(1, NumberOfBands+1):
		BandData = dataset.GetRasterBand(band).ReadAsArray()
		Prior[:,:,band-1] = BandData[InitRow:EndRow+1,0:cols]

	return Prior


def GetFileList(strInitDoY, DataDir):
	# Standard deviation and mean of the normal distribution to create the weighting factors	
	SD = 6
	mu = 0
	f = NormalDistribution(mu, SD)

	TimeWindow = 32
	Wing = TimeWindow / 2
	FileList = []
	InitDoY = int(strInitDoY)

	# Since the central DoY will be in the middle of the 16-day time period
	InitDoY = InitDoY + 8

	for Day in range(InitDoY - Wing, InitDoY + Wing):
		strDay = str(Day)
		if len(strDay)==1:
			strDay='00' + strDay
		elif len(strDay)==2:
			strDay='0' + strDay

		File = glob.glob(DataDir + '/AVH09C1.A????'+ strDay + '.N??.???.?????????????.tif')
	
		if len(File) == 1: FileList.append(File[0])

	FileList.sort()

	Year = numpy.zeros((len(FileList)), numpy.int16)
	DoY = numpy.zeros((len(FileList)), numpy.int16)
	Weigth = numpy.zeros((len(FileList)), numpy.float32)

	i = 0
	for File in FileList:
		# Get Year and DoY from filename
		YearOfObservation = os.path.basename(File).split('.')[1][1:5]
		DoYOfObservation = os.path.basename(File).split('.')[1][5:8]

		Year[i] = YearOfObservation
		DoY[i] = DoYOfObservation
		Weigth[i] = f(i - Wing) * (Wing - 1)

		i += 1

	return FileList, Year, DoY, Weigth


def GetDimensions(File):
	dataset = gdal.Open(File, GA_ReadOnly)
	# Usually the AVH09C1 layerstack should contain 10 bands
	rows, cols, NumberOfBands = dataset.RasterYSize, dataset.RasterXSize, dataset.RasterCount
	dataset = None

	return rows, cols, NumberOfBands

def BRDF_Inverter(Reflectance, Kernels, DoY, Weigth, NumberOfParameters, Prior):
	rows, cols, NumberOfBands, days = Reflectance.shape
	parameters = numpy.zeros((rows, cols, NumberOfBands, NumberOfParameters), numpy.float32)
	samples = numpy.zeros((rows, cols), numpy.float32)

	#(900, 7200, 10, 31)
	for j in range(0,cols):
		if j % 100 == 0:
			print "Processing columns", j
		for i in range(0,rows):
			# Extract temporal profile for broadbands vis, NIR  and short
			ReflectanceProfile = Reflectance[i,j,:,:]
			# Based on valid observations on vis
			NumberOfSamples = len(numpy.where( ReflectanceProfile[0] > 0)[0])

			if NumberOfSamples >= 7:
				IndexData = numpy.where(ReflectanceProfile[0] > 0)[0]
				w = Weigth[IndexData]
				refl = ReflectanceProfile[:,IndexData] * w

				K = numpy.ones([3, NumberOfSamples]) * w
				K[1,:] = Kernels[i,j,0,IndexData] * w
				K[2,:] = Kernels[i,j,1,IndexData] * w

				(P, rho_residuals, rank, svals) = lstsq( K.T, refl.T )
				parameters[i,j,:,:] = P
				samples[i,j] = w.sum()

				# Uncertainty for the f0 paramter is given by
				# :ref:`Lucht & Lewis`
				# s^2 =  MSE * M-1[0,0], where M^{-1} is
				# K*K^{T}
				M = numpy.dot ( K, K.T )
				unc_factor = numpy.linalg.inv ( M ) [0,0]
				unc = rho_residuals * unc_factor

				#if Prior[i,j][18] > 0:
				#	ipshell = embed()

	return parameters, samples


#--------------------------------------------------------------------------------#
from IPython import embed

strDoY = sys.argv[1]
DataDir = sys.argv[2]
PriorDataDir = "/Users/glopez/GCII/data/MODIS/MCD43C2/BroadBands"

FileList, Year, DoY, Weigth = GetFileList(strDoY, DataDir)

# From the first file get dimensions
rows, cols, NumberOfBands = GetDimensions(FileList[0])
# Overwrite number of bands with 3, only:
# vis, NIR and shortwave
NumberOfBands = 3
NumberOfParameters = 3

# Create aray where to store the composite
# all bands in LTDR plus DoY of observation, NumberOfSamples
parameters = numpy.zeros((rows, cols, NumberOfBands, NumberOfParameters), numpy.float32)
NumberOfSamples = numpy.zeros((rows, cols), numpy.float32)

# Depending on the processing system, the composite could be created storing ALL
# datasets in RAM, however for prototyping a tile-based processing will be implemented.
# 8 tiles will be the default setting.
NumberOfTiles = 8

for Tile in range(1,NumberOfTiles+1):
	InitRow = (Tile - 1) * (rows / NumberOfTiles)
	EndRow = (Tile * (rows / NumberOfTiles)) - 1
	print "Processing rows:", InitRow, "to", EndRow

	Prior = GetPrior(PriorDataDir, strDoY, InitRow, EndRow)

	# Create tmp layerstack
	NumberOfFiles = len(FileList)
	Reflectance = numpy.zeros(((rows / NumberOfTiles), cols, NumberOfBands, NumberOfFiles), numpy.float32)
	Kernels = numpy.zeros(((rows / NumberOfTiles), cols, 2, NumberOfFiles), numpy.float32)

	FileNumber = 1
	for File in FileList:
		print File
		dataset = gdal.Open(File, GA_ReadOnly)
		for band in range(1, NumberOfBands+1):
			BandData = dataset.GetRasterBand(band).ReadAsArray()
			if NumberOfTiles > 1:
				Reflectance[:,:,band-1,FileNumber-1] = BandData[InitRow:EndRow+1,0:cols]
			else:
				Reflectance[:,:,band-1,FileNumber-1] = BandData[:,:]
			BandData = None

		dataset = None

		KernelsFile = File[0:len(File)-3] + "kernels.tif"
		for band in range(1,3):
			dataset = gdal.Open(KernelsFile, GA_ReadOnly)
			BandData = dataset.GetRasterBand(band).ReadAsArray()
			if NumberOfTiles > 1:
				Kernels[:,:,band-1,FileNumber-1] = BandData[InitRow:EndRow+1,0:cols]
			else:
				Kernels[:,:,band-1,FileNumber-1] = BandData[:,:]

		BandData = None
		dataset = None

		FileNumber += 1

	print "Performing BRDF model inversion..."
	#ipshell = embed()

	tmpParameters, tmpNumberOfSamples = BRDF_Inverter(Reflectance, Kernels, DoY, Weigth, NumberOfParameters, Prior)
	parameters[InitRow:EndRow+1,0:cols,:,:] = tmpParameters
	NumberOfSamples[InitRow:EndRow+1,0:cols] = tmpNumberOfSamples
	tmpParameters = None
	tmpNumberOfSamples = None

print "Writing results to a file..."
format = "GTiff"
driver = gdal.GetDriverByName(format)
new_dataset = driver.Create( 'BRDF_Parameters.'+ strDoY +'.tif', cols, rows, (NumberOfBands*NumberOfParameters) + 1, GDT_Float32, ['COMPRESS=PACKBITS'])

k = 1
for i in range(1,NumberOfBands+1):
	for j in range(1,NumberOfParameters+1):
		new_dataset.GetRasterBand(k).WriteArray(parameters[:,:,j-1,i-1])
		k += 1

new_dataset.GetRasterBand((NumberOfBands*NumberOfParameters) + 1).WriteArray(NumberOfSamples[:,:])
new_dataset = None

