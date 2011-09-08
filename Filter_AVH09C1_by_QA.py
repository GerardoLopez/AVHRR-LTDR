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



