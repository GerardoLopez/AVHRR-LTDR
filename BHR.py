#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import sys

import osgeo.gdal as gdal
from osgeo.gdalconst import *
import numpy

BRDF = sys.argv[1]

#Get raster size
dataset = gdal.Open( BRDF, GA_ReadOnly )
ymax, xmax = dataset.RasterYSize, dataset.RasterXSize

NumberOfBands = 3
NumberOfParameters = 3

BHR = numpy.zeros((ymax, xmax, NumberOfBands), numpy.float32 )

i=0
#Get BRDF parameters of each band
for band_number in range (1, NumberOfBands*NumberOfParameters, NumberOfBands):
	f0 = dataset.GetRasterBand(band_number).ReadAsArray()
	f1 = dataset.GetRasterBand(band_number+1).ReadAsArray()
	f2 = dataset.GetRasterBand(band_number+2).ReadAsArray()

	BHR[:,:,i] = f0 + f1 * (0.189184) + f2 * (-1.377622)
	BHR[:,:,i] = numpy.where( ((BHR[:,:,i] < 0.0) | (BHR[:,:,i] > 1.0)), 0.0, BHR[:,:,i])

	i = i + 1


format = "ENVI"
driver = gdal.GetDriverByName(format)

new_dataset = driver.Create( sys.argv[1] + '.BHR.bin', xmax, ymax, NumberOfBands, GDT_Float32 )

for band_number in range (1, NumberOfBands+1):
	new_dataset.GetRasterBand(band_number).WriteArray(BHR[:,:,band_number-1])

new_dataset = None


