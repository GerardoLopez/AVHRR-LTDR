#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import sys

import osgeo.gdal as gdal
from osgeo.gdalconst import *
import numpy

# Read layerstack
print "Opening layerstack..."
dataset = gdal.Open( 'SDS_layerstack.img', GA_ReadOnly )
# Usually the AVH09C1 layerstack should contain 10 bands
rows = dataset.RasterYSize
cols = dataset.RasterXSize
NumberOfBands = dataset.RasterCount

layerstack = numpy.zeros((rows, cols, NumberOfBands), numpy.int16)
NumberOfBroadbands = 3
BB = numpy.zeros((rows, cols, NumberOfBroadbands), numpy.float32)
BB_error = numpy.zeros((rows, cols, NumberOfBroadbands), numpy.float32)

print "Reading bands..."
for i in range(NumberOfBands):
	layerstack[:,:,i] = dataset.GetRasterBand(i+1).ReadAsArray()

dataset = None

# Store angular information to calculate Kernels (class receives vza, sza, raa)
#   SZEN band 7
#   VZEN band 8
#   RELAZ band 9
VZA = layerstack[:,:,7] / 100.0
SZA = layerstack[:,:,6] / 100.0
RAA = layerstack[:,:,8] / 100.0

# QA bit-flag band is stored in band number 10 
QA = layerstack[:,:,9]
QA_binary = numpy.zeros((rows, cols), numpy.int16)
QA_error = numpy.zeros((rows, cols), numpy.float16)
QA_flags = [128, 160, 8320, 16512, -32640, -32608, -16256]
AbsoluteErrorBasedOnQA = [0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
#QA_flags = [128, 16512]

print "Screening data using QA flags..."
i=0
for flag in QA_flags:
	QA_binary = numpy.where(QA == flag, 1, QA_binary)
	QA_error = numpy.where(QA == flag, AbsoluteErrorBasedOnQA[i], QA_error)
	i += 1

for i in range(NumberOfBands):
	if i >= 0 and i <= 2:
		# Filter values with artifacts in reflectance bands
		layerstack[:,:,i] = numpy.where( (layerstack[:,:,i] < 0.0) | (layerstack[:,:,i] >= 10000), 0, layerstack[:,:,i])
	layerstack[:,:,i] = layerstack[:,:,i] * QA_binary

#===============================
# Narrow to broadband conversion
#===============================
# Conversion formulae and Residual Standard Error (RSE) taken from Liang, S. (2000)
# Narrowband to broadband conversions of land surface albedo I Algorithms
# Remote Sensing of Environment 76 (2000) 213-238

# Narrow bands
b1 = layerstack[:,:,0] / 10000.
b2 = layerstack[:,:,1] / 10000.

#---------------------Visible--------------------#
RSE_BB_visible = 0.0216
# BB_visible = 0.0074 + 0.5975*b1 + 0.4410*b1^2
BB[:,:,0] = 0.0074 + \
            (0.5975 * b1) + \
            (0.4410 * numpy.power(b1, 2))
BB[:,:,0] = numpy.where( (BB[:,:,0] < 0.0) | (BB[:,:,0] > 1.0), 0.0, BB[:,:,0] )
BB_error[:,:,0] = numpy.where( BB[:,:,0] > 0.0, (BB[:,:,0] * QA_error) + (BB[:,:,0] * RSE_BB_visible), BB_error[:,:,0] )

#---------------------NIR--------------------#
RSE_BB_NIR = 0.030
# BB_NIR = -1.4759*b1^2 - 0.6536*b2^2 + 1.8591*b1*b2 + 1.063*b2
BB[:,:,1] = (-1.4759 * numpy.power(b1, 2)) + \
            (-0.6536 * numpy.power(b2, 2)) + \
            (1.8591 * b1 * b2) + \
            (1.063 * b2)
BB[:,:,1] = numpy.where( (BB[:,:,1] < 0.0) | (BB[:,:,1] > 1.0), 0.0, BB[:,:,1] )
BB_error[:,:,1] = numpy.where( BB[:,:,1] > 0.0, (BB[:,:,1] * QA_error) + (BB[:,:,1] * RSE_BB_NIR), BB_error[:,:,1] )

#---------------------Shortwave--------------------#
RSE_BB_short = RSE_BB_visible + RSE_BB_NIR
# BB_short = -0.337*b1^2 - 0.270*b2^2 + 0.7074*b1*b2 + 0.2915*b1 + 0.5256*b2 + 0.0035
BB[:,:,2] = (-0.337 * numpy.power(b1, 2)) + \
            (-0.270 * numpy.power(b2, 2)) + \
            (0.7074 * b1 * b2) + \
            (0.2915 * b1) + \
            (0.5256 * b2) + \
            0.0035
BB[:,:,2] = numpy.where( (BB[:,:,2] < 0.0) | (BB[:,:,2] > 1.0), 0.0, BB[:,:,2] )
BB_error[:,:,2] = numpy.where( BB[:,:,2] > 0.0, (BB[:,:,2] * QA_error) + (BB[:,:,2] * RSE_BB_short), BB_error[:,:,2] )

# Save results to a ENVI binary file
format = "GTiff"
driver = gdal.GetDriverByName(format)
# Save to output file, 3 broadbands plus correspondent error, plus SREFL_CH3
NumberOfOutputBands = NumberOfBroadbands + NumberOfBroadbands + 1
new_dataset = driver.Create( 'SDS_layerstack_masked.tif', cols, rows, NumberOfOutputBands, GDT_Float32 , ['COMPRESS=PACKBITS'])

print "Saving data..."
BandNames=['BB_visible', 'BB_NIR', 'BB_short', 'BB_visible_error', 'BB_NIR_error', 'BB_short_error', 'SREFL_CH3']

for i in range(NumberOfBroadbands):
    new_dataset.GetRasterBand(i+1).WriteArray(BB[:,:,i] * QA_binary)
    new_dataset.SetMetadataItem("Band_" + str(i+1), BandNames[i])

for i in range(NumberOfBroadbands):
    new_dataset.GetRasterBand(NumberOfBroadbands+i+1).WriteArray(BB_error[:,:,i] * QA_binary)
    new_dataset.SetMetadataItem("Band_" + str(NumberOfBroadbands+i+1), BandNames[NumberOfBroadbands+i])

# Store SREFL_CH3 data
new_dataset.GetRasterBand(NumberOfOutputBands).WriteArray(layerstack[:,:,2])
new_dataset.SetMetadataItem("Band_" + str(NumberOfOutputBands), BandNames[NumberOfOutputBands-1])

new_dataset = None

# Let's tyde up a bit before calculating the kernels
BB = b1 = b2 = layerstack = None
layerstack = None
QA = QA_binary = None

#==================
# Calculate Kernels
#==================
print 'Calculating Ross-Thick and Li-Sparse kernels...'
import sys
sys.path.append('/Users/glopez/GCII/src/AVHRR-LTDR')
from Kernels import Kernels

# Kernels receives vza, sza, raa
kk = Kernels(VZA, SZA, RAA, \
             RossHS=False, RecipFlag=True, MODISSPARSE=True, \
             normalise=1, doIntegrals=False, LiType='Sparse', RossType='Thick')

Ross_Thick = numpy.zeros((rows,cols), numpy.float32 )
Li_Sparse = numpy.zeros((rows, cols), numpy.float32 )

Ross_Thick = kk.Ross.reshape((rows,cols))
Li_Sparse = kk.Li.reshape((rows,cols))

kk = None

# Save kernels in a different file to keep all 32bit floating point information
new_dataset = driver.Create( 'SDS_layerstack_kernels_masked.tif', cols, rows, 2, GDT_Float32, ['COMPRESS=PACKBITS'] )
new_dataset.GetRasterBand(1).WriteArray(Ross_Thick)
new_dataset.GetRasterBand(2).WriteArray(Li_Sparse)
new_dataset = None

#from IPython import embed
#ipshell = embed()

