#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import sys

import osgeo.gdal as gdal
from osgeo.gdalconst import *
import numpy

BHR = sys.argv[1]

#Get raster size
dataset = gdal.Open( BHR, GA_ReadOnly )
ymax, xmax = dataset.RasterYSize, dataset.RasterXSize

NumberOfBands = 3

NDVI = numpy.zeros((ymax, xmax), numpy.float32 )

vis = dataset.GetRasterBand(1).ReadAsArray()
nir = dataset.GetRasterBand(2).ReadAsArray()

NDVI = numpy.divide((nir - vis), (nir + vis))
NDVI = numpy.where(numpy.isnan(NDVI)==True, 0, NDVI)

format = "ENVI"
driver = gdal.GetDriverByName(format)

new_dataset = driver.Create( sys.argv[1] + '.NDVI.bin', xmax, ymax, 1, GDT_Float32 )
new_dataset.GetRasterBand(1).WriteArray(NDVI[:,:])
new_dataset = None


