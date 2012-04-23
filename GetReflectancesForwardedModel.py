#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

from BRDF import *
import sys

def GetModelParameters(BRDF_ParametersFile):
	#Get raster size
	rows, cols, NumberOfBands = GetDimensions(BRDF_ParametersFile)

	NumberOfBands = 3
	NumberOfParameters = 3

	Parameters = numpy.zeros((rows,cols,NumberOfBands * NumberOfParameters), numpy.float32)

	#Get BRDF parameters of each wave band
	dataset = gdal.Open(BRDF_ParametersFile, GA_ReadOnly)
	for band_number in range (NumberOfBands*NumberOfParameters):
		Parameters[:,:,band_number] = dataset.GetRasterBand(band_number+1).ReadAsArray()

	NSamples = dataset.GetRasterBand(19).ReadAsArray()

	dataset = None

	return Parameters, NSamples

def GetGlobalKernels(KernelsFile):
	#Get raster size
	rows, cols, NumberOfBands = GetDimensions(KernelsFile)

	# Get the Geometric and Volumetric kernels
	Kernels = numpy.zeros((rows,cols,2), numpy.float32)

	dataset = gdal.Open(KernelsFile, GA_ReadOnly)
	Kernels[:,:,0] = dataset.GetRasterBand(1).ReadAsArray()
	Kernels[:,:,1] = dataset.GetRasterBand(2).ReadAsArray()

	return Kernels


def GetReflectancesForwardedModel(Kernels, Parameters, Mask, NSamples):
	'''
	Get predicted reflectances based on model:
	refl = f0 + f1*Kvol + f2*KGeo
	Where:
		f0 = Isotropic parameter
		f1 = Volumetric parameter
		f2 = Geometric parameter

		KVol = Ross-Thick (volumetric kernel)
		KGeo = Li-Sparse (geometric kernel)

	Uncertainty in predicted reflectance
		Ra = Predicted reflectances from the model           PredictedReflectances
		invCa = the uncertainty inverse matrix in the model  Parameters.MData

		Robs = the observed reflectance                      Reflectance
		Cobs = the uncertainty matrix in reflectances        Cinv

		Xa^2 = (Robs - Ra)^T (Cobs + kernels^T invCa kernels) (Robs - Ra)

	'''

	nWaveBands = 3

	# Get predicted reflectances
	VIS = Parameters[:,:,0] + (Parameters[:,:,1] * Kernels[:,:,0] ) + (Parameters[:,:,2] * Kernels[:,:,1])
	NIR = Parameters[:,:,3] + (Parameters[:,:,4] * Kernels[:,:,0] ) + (Parameters[:,:,5] * Kernels[:,:,1])
	SW =  Parameters[:,:,6] + (Parameters[:,:,7] * Kernels[:,:,0] ) + (Parameters[:,:,8] * Kernels[:,:,1])

	PredictedReflectances = numpy.zeros((nWaveBands, Parameters.shape[0], Parameters.shape[1]), numpy.float32)
	#Mask data
	PredictedReflectances[0,:,:] = numpy.where((Mask == 0) | (NSamples==0), 0.0, VIS)
	PredictedReflectances[1,:,:] = numpy.where((Mask == 0) | (NSamples==0), 0.0, NIR)
	PredictedReflectances[2,:,:] = numpy.where((Mask == 0) | (NSamples==0), 0.0, SW)

    #Xa = (Robs - Ra)^T (Cobs + kernels^T invCa kernels) (Robs - Ra)

    ##PredictedReflectances_SD = numpy.zeros((nWaveBands, parameters.shape[1], parameters.shape[2]), numpy.float32)
    ##Xa = numpy.zeros((parameters.shape[1], parameters.shape[2]), numpy.float32)
    ##Robs = BBDR
    ##Ra = PredictedReflectances
    ##Cobs = Cinv
    ##invCa = Prior.MData

	##for columns in range(0,parameters.shape[1]):
    ##    for rows in range(0, parameters.shape[2]):
    ##        # Only examine pixels where observed reflectances are greater than predicted reflectances in 3 BB
    ##        if Mask[columns,rows] <> 0 and numpy.where((Robs[:,columns,rows] - Ra[:,columns,rows]>0))[0].shape[0]==3 :
                #Xa = (Robs - Ra)^T (Cobs + kernels^T invCa kernels) (Robs - Ra)
    ##            k = numpy.zeros((9,3), numpy.float32)
    ##            k[0:3,0] = kernels[0:3,columns,rows]
    ##            k[3:6,1] = kernels[3:6,columns,rows]
    ##            k[6:9,2] = kernels[6:9,columns,rows]

    ##            Xa[columns,rows] = numpy.matrix(Robs[:,columns,rows] - Ra[:,columns,rows]) * \
    ##                          (numpy.matrix(Cobs[:,:,columns,rows]) + numpy.matrix(k).T * numpy.matrix(invCa[:,:,columns,rows]) * numpy.matrix(k)) * \
    ##                           numpy.matrix(Robs[:,columns,rows] - Ra[:,columns,rows]).T

    ##Xa = numpy.sqrt(Xa)

	#return ReturnReflectancesForwardedModel(PredictedReflectances)
	return PredictedReflectances

class ReturnReflectancesForwardedModel(object):
	def __init__(self, PredictedReflectances):
		self.PredictedReflectances = PredictedReflectances
		#self.Xa = Xa

# -------------------------------------------------------------------------------- #
from IPython import embed

BRDF_ParametersFile = sys.argv[1]
ReflectancesFile = sys.argv[2]
KernelsFile = sys.argv[3]

Parameters, NSamples = GetModelParameters(BRDF_ParametersFile)
InitRow = 0
EndRow = Parameters.shape[0] - 1

Reflectance = GetOnlyReflectances(ReflectancesFile, InitRow, EndRow)
Mask = numpy.where(Reflectance[:,:,0] == 0.0, 0, 1)
Kernels = GetGlobalKernels(KernelsFile)

ModeledReflectances = GetReflectancesForwardedModel(Kernels, Parameters, Mask, NSamples)

format = "GTiff"
driver = gdal.GetDriverByName(format)

new_dataset = driver.Create( 'ModeledReflectances.tif', Parameters.shape[1], Parameters.shape[0], 3, GDT_Float32, ['COMPRESS=PACKBITS'] )

new_dataset.GetRasterBand(1).WriteArray(ModeledReflectances[0,:,:])
new_dataset.GetRasterBand(2).WriteArray(ModeledReflectances[1,:,:])
new_dataset.GetRasterBand(3).WriteArray(ModeledReflectances[2,:,:])
new_dataset = None

#ipshell = embed()


