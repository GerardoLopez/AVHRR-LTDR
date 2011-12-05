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

print "Reading bands..."
for i in range(1,NumberOfBands+1):
	layerstack[:,:,i-1] = dataset.GetRasterBand(i).ReadAsArray()

dataset = None

# QA bit-flag band is stored in band number 10 
QA = layerstack[:,:,9]
QA_binary = numpy.zeros((rows, cols), numpy.int16)
#QA_flags = [128, 16512, 8320]
QA_flags = [128, 16512]

print "Screening data using QA flags..."
for flag in QA_flags:
	QA_binary = numpy.where(QA == flag, 1, QA_binary)

for i in range(1,NumberOfBands+1):
	layerstack[:,:,i-1] = layerstack[:,:,i-1] * QA_binary

# Calculate Kernels
#   SZEN band 7
#   VZEN band 8
#   RELAZ band 9
print 'Calculating Ross-Thick and Li-Sparse kernels...'
import sys
sys.path.append('/Users/glopez/GCII/src/AVHRR-LTDR')
from Kernels import *

# Kernels receives vza, sza, raa
kk = Kernels(layerstack[:,:,7]/100.0 ,layerstack[:,:,6]/100.0 , layerstack[:,:,8]/100.0, RossHS=False, MODISSPARSE=True, RecipFlag=True, normalise=1, doIntegrals=False, LiType='Sparse', RossType='Thick')

Ross_Thick = numpy.zeros((rows,cols), numpy.float32 )
Li_Sparse = numpy.zeros((rows, cols), numpy.float32 )

Ross_Thick = kk.Ross.reshape((rows,cols))
Li_Sparse = kk.Li.reshape((rows,cols))

kk = None

# Save results to a ENVI binary file
format = "GTiff"
driver = gdal.GetDriverByName(format)
new_dataset = driver.Create( 'SDS_layerstack_masked.tif', cols, rows, NumberOfBands, GDT_Int16 , ['COMPRESS=PACKBITS'])

print "Saving data..."
BandNames=['SREFL_CH1', 'SREFL_CH2', 'SREFL_CH3', 'BT_CH3', 'BT_CH4', 'BT_CH5', 'SZEN', 'VZEN', 'RELAZ', 'QA' ]

for i in range(1,NumberOfBands+1):
	new_dataset.GetRasterBand(i).WriteArray(layerstack[:,:,i-1])
	new_dataset.SetMetadataItem("Band_" + str(i), BandNames[i-1])

new_dataset = None

# Save kernels in a different file to keep all 32bit floating point information
new_dataset = driver.Create( 'SDS_layerstack_kernels_masked.tif', cols, rows, 2, GDT_Float32, ['COMPRESS=PACKBITS'] )
new_dataset.GetRasterBand(1).WriteArray(Ross_Thick)
new_dataset.GetRasterBand(2).WriteArray(Li_Sparse)
new_dataset = None

#from IPython import embed
#ipshell = embed()

