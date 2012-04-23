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

def GetPixel(File, Year, DoY, xmin=1, ymin=1, xmax=1, ymax=1):
    import sys
    import glob

    try:
        filename = glob.glob(File)
		# Fpar_1km Lai_1km FparLai_QC FparExtra_QC FparStdDev_1km LaiStdDev_1km
        dataset = gdal.Open( filename[0], GA_ReadOnly )
    except:
        print "Error:", sys.exc_info()[0]
        exit(-1)

    if xmin == ymin == xmax == ymax == 1:
        Xmin = Ymin = 0
        Xmax = dataset.GetSubDatasets()[0][1].split()[0].split('x')[0].split('[')[1]
        Ymax = dataset.GetSubDatasets()[0][1].split()[0].split('x')[1].split(']')[0]
        BandCount = len(dataset.GetSubDatasets())
        xsize = (Xmax-Xmin)
        ysize = (Ymax-Ymin)

    else:
        BandCount = len(dataset.GetSubDatasets())
        xsize = (xmax-xmin) + 1
        ysize = (ymax-ymin) + 1
        Xmin = xmin - 1
        Ymin = ymin - 1

    SubDatasets = dataset.GetSubDatasets()
    dataset = None

    #print "Opening file", filename[0]

    # Load data
	# Change data type according to input data, e.g., float32
    data = numpy.zeros((xsize, ysize, BandCount), numpy.int16)

    for i in range(BandCount):
        dataset = gdal.Open( SubDatasets[i][0], GA_ReadOnly )
        data[:,:,i] = dataset.GetRasterBand(1).ReadAsArray(Xmin, Ymin, xsize, ysize)
        dataset = None

    strData = str(int(Year)) + ',' + str(int(DoY))
    for band in range(BandCount):
        # If data type is float change int to flot
        strData = strData + ',' + str(int(data[:,:,band]))

    print strData


import sys

File = sys.argv[1]

xmin = int(sys.argv[2])
ymin = int(sys.argv[3])
xmax = int(sys.argv[4])
ymax = int(sys.argv[5])

Year = sys.argv[6]
DoY = sys.argv[7]

results = GetPixel(File, Year, DoY, xmin, ymin, xmax, ymax)
