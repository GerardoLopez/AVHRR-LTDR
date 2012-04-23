#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import glob
import os
import sys
import time

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

def GetOnlyReflectances(ReflectancesFile, InitRow, EndRow):
	rows, cols, NumberOfBands = GetDimensions(ReflectancesFile)
	rows = (EndRow-InitRow)+1

	# Band 1 = vis    Band 2 = NIR    Band 3 = shortwave
	NumberOfBands = 3

	Reflectance = numpy.zeros((rows, cols, NumberOfBands), numpy.float32)

	dataset = gdal.Open(ReflectancesFile, GA_ReadOnly)
	for i in range(NumberOfBands):
		BandData = dataset.GetRasterBand(i+1).ReadAsArray(0,InitRow, cols ,rows)
		Reflectance[:,:,i] = BandData

	dataset = None

	return Reflectance


def GetReflectances(ReflectancesFile, KernelsFile, Weigth, InitRow, EndRow):
	rows, cols, NumberOfBands = GetDimensions(KernelsFile)
	rows = (EndRow-InitRow)+1

	# Band 1 = vis    Band 2 = NIR    Band 3 = shortwave
	# Band 4 = vis SD Band 5 = NIR SD Band 6 = shortwave SD
	NumberOfBands = 3

	# BBDR matrix size -> rows x columns x NumberOfBands x 1
	Reflectance = numpy.zeros((rows, cols, NumberOfBands, 1), numpy.float32)
	ReflectanceSD = numpy.zeros((rows, cols, NumberOfBands), numpy.float32)
	NumberOfSamples = numpy.zeros((rows, cols), numpy.float32)

	dataset = gdal.Open(ReflectancesFile, GA_ReadOnly)
	for i in range(NumberOfBands):

		#BandData = dataset.GetRasterBand(i+1).ReadAsArray()
		#Reflectance[:,:,i,0] = BandData[InitRow:EndRow+1,0:cols]
		BandData = dataset.GetRasterBand(i+1).ReadAsArray(0,InitRow, cols ,rows)
		Reflectance[:,:,i,0] = BandData

		#BandData = dataset.GetRasterBand(i+NumberOfBands+1).ReadAsArray()
		#ReflectanceSD[:,:,i] = BandData[InitRow:EndRow+1,0:cols]
		BandData = dataset.GetRasterBand(i+NumberOfBands+1).ReadAsArray(0,InitRow, cols ,rows)
		ReflectanceSD[:,:,i] = BandData

	dataset = None

	#-------------------------------------
	# Build C -- coveriance matrix for obs
	#-------------------------------------
	# C in a symetric matrix form of NumberOfBands * NumberOfBands
	C = numpy.zeros((rows, cols, NumberOfBands, NumberOfBands), numpy.float32)
	Cinv = numpy.zeros((rows, cols, NumberOfBands, NumberOfBands), numpy.float32)

	for i in range(NumberOfBands):
		C[:,:,i,i] = ReflectanceSD[:,:,i] * ReflectanceSD[:,:,i]

	# Create matrices: M, V and E
	# M = Kernels^T C^-1 Kernels
	# V = Kernels^T C^-1 Reflectance
	# E = Reflectance^T C^-1 Reflectance
	M = numpy.zeros((rows, cols, NumberOfBands*NumberOfBands, NumberOfBands*NumberOfBands), numpy.float32)
	V = numpy.zeros((rows, cols, NumberOfBands*NumberOfBands), numpy.float32)
	E = numpy.zeros((rows, cols), numpy.float32)

	# Get Kernels
	Kernels = GetKernels(KernelsFile, InitRow, EndRow)

	for j in range(0,cols):
		for i in range(0,rows):
			if numpy.where( (Reflectance[i,j,:]>0.0) & (Reflectance[i,j,:]<=1.0) )[0].shape[0] == NumberOfBands and \
               numpy.where( (ReflectanceSD[i,j,:]>=0.0) & (ReflectanceSD[i,j,:]<=1.0) )[0].shape[0] == NumberOfBands:

				Cinv[i,j,:,:] = numpy.matrix(C[i,j,:,:]).I
				M[i,j,:,:] = numpy.matrix(Kernels[i,j,:,:]).T * numpy.matrix(Cinv[i,j,:,:]) * numpy.matrix(Kernels[i,j,:,:])
				# Multiply only using lead diagonal of Cinv, additionally transpose the result to store V as a 1 x 9 vector
				V[i,j,:] = (numpy.matrix(Kernels[i,j,:,:]).T * numpy.diagflat(numpy.diagonal(Cinv[i,j,:,:])) * Reflectance[i,j,:,:]).T
				E[i,j] = numpy.matrix(Reflectance[i,j,:,:]).T * numpy.matrix(Cinv[i,j,:,:]) * Reflectance[i,j,:,:]
				NumberOfSamples[i,j] = Weigth

	return ReturnGetReflectances(Reflectance, M, V, E, NumberOfSamples, Kernels)


class ReturnGetReflectances(object):
    def __init__(self, Reflectance, M, V, E, NumberOfSamples, Kernels):
		self.Reflectance = Reflectance
		self.M = M
		self.V = V
		self.E = E
		self.NumberOfSamples = NumberOfSamples
		self.Kernels = Kernels


def GetKernels(KernelsFile, InitRow, EndRow):
	rows, cols, NumberOfBands = GetDimensions(KernelsFile)
	rows = (EndRow-InitRow)+1

	NumberOfBands = 3
	NumberOfKernels = NumberOfBands * 3

	# kernels matrix size -> rows x columns = 3 x 9
	Kernels = numpy.zeros((rows, cols, NumberOfBands, NumberOfKernels), numpy.float32)
    # Isotropic kernels = 1
	Kernels[:,:,0,0] = Kernels[:,:,1,3] = Kernels[:,:,2,6] = 1.0

	for i in range(1,NumberOfBands):
		dataset = gdal.Open(KernelsFile, GA_ReadOnly)
		BandData = dataset.GetRasterBand(i).ReadAsArray(0,InitRow, cols ,rows)
		Kernels[:,:,0,i] = Kernels[:,:,1,i+NumberOfBands] = Kernels[:,:,2,i+NumberOfBands*2] = BandData

	dataset = BandData = None

	return Kernels


def GetPrior(PriorDataDir, strDoY, InitRow, EndRow):
	PriorFile = PriorDataDir + '/MCD43C2.Prior.' + strDoY + '.img'
	print "Opening prior", PriorFile
	rows, cols, NumberOfBands = GetDimensions(PriorFile)

	# Overwright NumberOfBands to the number of wavebands
	NumberOfBands = NumberOfParameters = 3

	rows = (EndRow-InitRow)+1
	Prior = numpy.zeros((rows,cols,NumberOfBands * NumberOfParameters), numpy.float32)
	PriorVariance = numpy.zeros((rows,cols,NumberOfBands * NumberOfParameters), numpy.float32)
	Mask = numpy.zeros((rows,cols), numpy.int8)

	C = numpy.zeros((rows, cols, NumberOfBands * NumberOfParameters, NumberOfBands * NumberOfParameters), numpy.float32)
	Cinv = numpy.zeros((rows, cols, NumberOfBands * NumberOfParameters, NumberOfBands * NumberOfParameters), numpy.float32)
	CinvF = numpy.zeros((rows, cols, NumberOfBands * NumberOfParameters), numpy.float32) # Matrix to store C^-1 * Fpr

	dataset = gdal.Open(PriorFile, GA_ReadOnly)

	for i in range(NumberOfBands * NumberOfParameters):
		BandData = dataset.GetRasterBand(i+1).ReadAsArray(0,InitRow, cols ,rows)
		Prior[:,:,i] = BandData

		BandData = dataset.GetRasterBand(i + (NumberOfBands * NumberOfParameters) + 1).ReadAsArray(0,InitRow, cols ,rows)
		# Could be the case that the covariance is 0 or very small but there are samples, make variance = 1
		PriorVariance[:,:,i] = numpy.where(BandData[:,:] <= 1.0e-8, 1.0, BandData[:,:])
		C[:,:,i,i] = PriorVariance[:,:,i]

	BandData = dataset = None

	# Calculate C inverse
	for j in range(0,cols):
		for i in range(0,rows):
			# Check that al least the isotropic parameters have values
			if numpy.where( (Prior[i,j,[0,3,6]]>0.0) & (Prior[i,j,[0,3,6]]<=1.0) )[0].shape[0] == 3:

				try:
					Cinv[i,j,:,:] = numpy.matrix(C[i,j,:,:]).I
				except numpy.linalg.LinAlgError:
					indices = numpy.where(PriorVariance[i,j,:] == 0.0)[0]
					PriorVariance[i,j,indices] = 1.0
					C[i,j,indices,indices] = 1.0
					Cinv[i,j,:,:] = numpy.matrix(C[i,j,:,:]).I

				for k in range(NumberOfBands * NumberOfParameters):
					CinvF[i,j,k] = Cinv[i,j,k,k] * Prior[i,j,k]

				Mask[i,j] = 1

	return ReturnPriorData(Cinv, CinvF, Prior, PriorVariance, Mask)


class ReturnPriorData(object):
	def __init__(self, M, V, Parameters, ParametersVariance, Mask):
		self.M = M
		self.V = V
		self.Parameters = Parameters
		self.ParametersVariance = ParametersVariance
		self.Mask = Mask

def IsLeapYear(Year):
	if Year % 4 == 0:
		if Year % 100 == 0:
			if Year % 400 == 0:
				return True
			else:
				return False
		else:
			return True
	else:
		return False

def GetFileList(strInitDoY, DataDir, Year):
	# Standard deviation and mean of the normal distribution to create the weighting factors	
	SD = 6
	mu = 0
	f = NormalDistribution(mu, SD)

	TimeWindow = 32
	Wing = TimeWindow / 2
	FileList = []
	InitDoY = int(strInitDoY)

	# Since the central DoY will be in the middle of the 16-day time period
	InitDoY += 8

	for Day in range(InitDoY - Wing, InitDoY + Wing):

		# Some DoY could be negative
		# Some DoY could be greater than 365/366
		#     Assign the right DoY taking into account if Year is a leap-year
		if Day <= 0:
			if IsLeapYear(Year):
				Day = 366 + Day
			else:
				Day = 365 + Day
			strYear = str(Year - 1)

		elif Day >= 366 and IsLeapYear(Year):
			Day = Day - 365

		elif Day >= 367:
			if IsLeapYear(Year):
				Day = Day - 366
			else:
				Day = Day - 365
			strYear = str(Year + 1)
		else:
			strYear = str(Year)
			
		strDay = str(Day)
		if len(strDay)==1:
			strDay='00' + strDay
		elif len(strDay)==2:
			strDay='0' + strDay

		File = glob.glob(DataDir + '/' + strYear + '/AVH09C1.A'+ strYear + strDay + '.N??.???.?????????????.tif')
	
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

def BRDF_Inverter(M, V, E, Prior, NumberOfBands, NumberOfParameters, NumberOfSamples):
	rows, cols = E.shape
	Parameters = numpy.zeros((rows, cols, NumberOfBands, NumberOfParameters), numpy.float32)
	ParametersVariance = numpy.zeros((rows, cols, NumberOfBands, NumberOfParameters), numpy.float32)

	for j in range(0,cols):
		if j % 100 == 0:
			print "Processing columns", j
		for i in range(0,rows):
			if (Prior.Mask[i,j] == 1 and NumberOfSamples[i,j] > 0.0) or (NumberOfSamples[i,j] > 0.0):
				M_inversion = M[i,j,:,:] + Prior.M[i,j,:,:]
				V_inversion = V[i,j,:] + Prior.V[i,j,:]

				(P, rho_residuals, rank, svals) = lstsq(M_inversion, V_inversion)
				# TO DO: check under what conditions the BRDF params are negative
				# If negative, use Prior parameters
				indices = numpy.where(P < 0.0)
				P[indices] = P[indices] = Prior.Parameters[i,j,indices]

				Parameters[i,j,:,:] = P.reshape(NumberOfBands,NumberOfParameters)

				try:
					ParametersVariance[i,j,:,:] = numpy.diagonal(numpy.matrix(M_inversion).I).reshape(NumberOfBands,NumberOfParameters)
				except numpy.linalg.LinAlgError:
					# Set all off-diagonal elements of M_inversion to 0
					M_inversion_tmp = numpy.zeros((M_inversion.shape), numpy.float32)
					for k in range(0,9):
						M_inversion_tmp[k,k] = M_inversion[k,k]

					ParametersVariance[i,j,:,:] = numpy.diagonal(numpy.matrix(M_inversion_tmp).I).reshape(NumberOfBands,NumberOfParameters)
					M_inversion_tmp = None
					#ipshell = embed()
			else:
				# If there are no samples at all, use only prior data
				Parameters[i,j,:,:] = Prior.Parameters[i,j,:].reshape(NumberOfBands,NumberOfParameters)
				ParametersVariance[i,j,:,:] = Prior.ParametersVariance[i,j,:].reshape(NumberOfBands,NumberOfParameters)

	return Parameters, ParametersVariance

