#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

"""
Extract data for a single pixel

Authors:   Gerardo Lopez-Saldana <GerardoLopez@isa.utl.pt>
"""
import os

try:
  import numpy
except ImportError:
  print 'Numpy is not installed'
  exit(-1)

try:
  import osgeo.gdal as gdal
  from osgeo.gdalconst import *
except ImportError:
  print 'GDAL is not installed.'
  exit(-1)

def GetPixel(File, DoY, xmin=1, ymin=1, xmax=1, ymax=1):
    import sys
    import glob

    try:
        filename = glob.glob(File)
        dataset = gdal.Open( filename[0], GA_ReadOnly )
    except:
        print "Error:", sys.exc_info()[0]
        exit(-1)

    if xmin == ymin == xmax == ymax == 1:
        Xmin = Ymin = 0
        Xmax, Ymax, BandCount = dataset.RasterXSize, dataset.RasterYSize, dataset.RasterCount
        xsize = (Xmax-Xmin)
        ysize = (Ymax-Ymin)

    else:
        BandCount = dataset.RasterCount
        xsize = (xmax-xmin) + 1
        ysize = (ymax-ymin) + 1
        Xmin = xmin - 1
        Ymin = ymin - 1

    #print "Opening file", filename[0]

    # Load data
	# Change data type according to input data, e.g., float32
	data = numpy.zeros((xsize, ysize, BandCount), numpy.int16)

	for band in range(BandCount):
		data[:,:,band] = dataset.GetRasterBand(band+1).ReadAsArray(Xmin, Ymin, xsize, ysize)

	strData = str(float(DoY))
	for band in range(BandCount):
		# If data type is float change int to flot
		strData = strData + ',' + str(float(data[:,:,band]))

	print strData


import sys

File = sys.argv[1]

xmin = int(sys.argv[2])
ymin = int(sys.argv[3])
xmax = int(sys.argv[4])
ymax = int(sys.argv[5])

DoY = sys.argv[6]

results = GetPixel(File, DoY, xmin, ymin, xmax, ymax)
