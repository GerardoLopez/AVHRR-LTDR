#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import sys
import osgeo.gdal as gdal
from osgeo.gdalconst import *
import numpy

FireObsFilename = sys.argv[1]

InitDay = float(sys.argv[2])
strInitDay = sys.argv[2]
#e.g. 19980101

EndDay = float(sys.argv[3])
strEndDay = sys.argv[3]
#e.g. 19980116

# File with the master (ULX, ULY)
Master = '/Users/glopez/GCII/data/MODIS/MCD43C2/BroadBands/Master/Master.img'
dataset = gdal.Open( Master, GA_ReadOnly )
ymax, xmax = dataset.RasterYSize, dataset.RasterXSize
# FireClusters will have two bands, the fire count and the date of last fire observation
FireClusters = numpy.zeros((ymax, xmax, 2), numpy.int32 )
# Read GeoTransform
Projection = dataset.GetProjection()
GeoMatrix = dataset.GetGeoTransform()
(success, inv_geometrix) = gdal.InvGeoTransform(GeoMatrix)


# Open and read fire observations
# Longitude, Latitude, Month, Day, Julian date and Year
f = open(FireObsFilename)
fires = f.readlines()
f.close()

for fire in fires:
	FireData = fire.split()
	Year = FireData[5]
	if len(FireData[2]) == 1:
		Month = "0" + FireData[2]
	else:
		Month = FireData[2]

	if len(FireData[3]) == 1:
		Day = "0" + FireData[3]
	else:
		Day = FireData[3]

	date = Year + Month + Day
	latitude = FireData[0]
	longitude = FireData[1]

	if float(date) >= InitDay and float(date) <= EndDay:
		X = float(longitude)
		Y = float(latitude)

		# Calculate corresponding pixel coordinates (row/column) for a specific Lat/Lon location
		x = int(inv_geometrix[0] + inv_geometrix[1] * X + inv_geometrix[2] * Y)
		y = int(inv_geometrix[3] + inv_geometrix[4] * X + inv_geometrix[5] * Y)
		#print x, y, date
		FireClusters[y+1, x+1, 0] = FireClusters[y+1, x+1, 0] + 1
		FireClusters[y+1, x+1, 1] = long(date)

filename = strInitDay + '_' + strEndDay + '.txt'
f = open(filename, 'w')

# Find pixels where there are more than 4 ATSR-derived fire observations
indices = numpy.where(FireClusters[:,:,0] >= 4)
for i in range(len(indices[0])):
	date = str(int(FireClusters[indices[0][i], indices[1][i], 1]))
	year = date[0:4]
	month = date[4:6]
	day = date[6:8]
	record = str(indices[1][i]+1) + " " + str(indices[0][i]+1) + " " + str(int(FireClusters[indices[0][i], indices[1][i], 0])) + " "+ day + " " + month + " " + year + "\n"
	f.write(record)

f.close()

#from IPython import embed
#ipshell = embed()



# Save file
#print "Writing results to a file..."
#format = "GTiff"
#driver = gdal.GetDriverByName(format)
#new_dataset = driver.Create( strInitDay + '_' + strEndDay +'.tif', xmax, ymax, 2, GDT_Float32)
#new_dataset.SetGeoTransform(GeoMatrix)
#new_dataset.SetProjection(Projection)
#new_dataset.GetRasterBand(1).WriteArray(FireClusters[:,:,0])
#new_dataset.GetRasterBand(2).WriteArray(FireClusters[:,:,1])
#new_dataset = None



