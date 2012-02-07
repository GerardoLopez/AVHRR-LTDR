#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import BRDF.py

def GetModelParameters(BRDF_ParametersFile)
	#Get raster size
	rows, cols, NumberOfBands = GetDimensions(BRDF_ParametersFile)

	NumberOfBands = 3
	NumberOfParameters = 3

	Parameters = numpy.zeros((rows,cols,NumberOfBands * NumberOfParameters), numpy.float32)

	#Get BRDF parameters of each wave band
	dataset = gdal.Open(BRDF_ParametersFile, GA_ReadOnly)
	for band_number in range (NumberOfBands*NumberOfParameters):
		Parameters[:,:,i] = dataset.GetRasterBand(i+1).ReadAsArray()

	return Parameters


def GetReflectancesForwardedModel(Kernels, Parameters):
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

    #Mask data
    Parameters = numpy.where(Kernels==0, 0.0, Parameters)
    parameters_SD = numpy.where(kernels==0, 0.0, Prior.Parameters_SD)

    # Get predicted reflectances
    VIS = parameters[0,:,:] + (parameters[1,:,:] * kernels[1,:,:] ) + (parameters[2,:,:] * kernels[2,:,:])
    NIR = parameters[3,:,:] + (parameters[4,:,:] * kernels[4,:,:] ) + (parameters[5,:,:] * kernels[5,:,:])
    SW =  parameters[6,:,:] + (parameters[7,:,:] * kernels[7,:,:] ) + (parameters[8,:,:] * kernels[8,:,:])

    PredictedReflectances = numpy.zeros((nWaveBands, parameters.shape[1], parameters.shape[2]), numpy.float32)
    PredictedReflectances[0,:,:] = VIS
    PredictedReflectances[1,:,:] = NIR
    PredictedReflectances[2,:,:] = SW

    #Xa = (Robs - Ra)^T (Cobs + kernels^T invCa kernels) (Robs - Ra)
    PredictedReflectances_SD = numpy.zeros((nWaveBands, parameters.shape[1], parameters.shape[2]), numpy.float32)
    Xa = numpy.zeros((parameters.shape[1], parameters.shape[2]), numpy.float32)
    Robs = BBDR
    Ra = PredictedReflectances
    Cobs = Cinv
    invCa = Prior.MData

for columns in range(0,parameters.shape[1]):
        for rows in range(0, parameters.shape[2]):
            # Only examine pixels where observed reflectances are greater than predicted reflectances in 3 BB
            if Mask[columns,rows] <> 0 and numpy.where((Robs[:,columns,rows] - Ra[:,columns,rows]>0))[0].shape[0]==3 :
                #Xa = (Robs - Ra)^T (Cobs + kernels^T invCa kernels) (Robs - Ra)
                k = numpy.zeros((9,3), numpy.float32)
                k[0:3,0] = kernels[0:3,columns,rows]
                k[3:6,1] = kernels[3:6,columns,rows]
                k[6:9,2] = kernels[6:9,columns,rows]

                Xa[columns,rows] = numpy.matrix(Robs[:,columns,rows] - Ra[:,columns,rows]) * \
                              (numpy.matrix(Cobs[:,:,columns,rows]) + numpy.matrix(k).T * numpy.matrix(invCa[:,:,columns,rows]) * numpy.matrix(k)) * \
                               numpy.matrix(Robs[:,columns,rows] - Ra[:,columns,rows]).T

    Xa = numpy.sqrt(Xa)

    return ReturnReflectancesForwardedModel(PredictedReflectances, Xa)

class ReturnReflectancesForwardedModel(object):
    def __init__(self, PredictedReflectances, Xa):
        self.PredictedReflectances = PredictedReflectances
        self.Xa = Xa
