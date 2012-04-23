#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import sys
import shapefile


filename = sys.argv[1]
InitDay = float(sys.argv[2])
strInitDay = sys.argv[2]
EndDay = float(sys.argv[3])
strEndDay = sys.argv[3]


f = open(filename)
fires = f.readlines()

w = shapefile.Writer(shapefile.POINT)
w.field('Date')
w.field('Orbit')
w.field('Time')
w.field('Latitude')
w.field('Longitude')

for fire in fires:
	FireData = fire.split()
	date = FireData[0]
	orbit = FireData[1]
	time = FireData[2]
	latitude = FireData[3]
	longitude = FireData[4]

	if float(date) >= InitDay and float(date) <= EndDay:
		print fire
		w.point(float(longitude), float(latitude))
		w.record(long(date), long(orbit), float(time), float(latitude), float(longitude))


w.save( strInitDay + '_' + strEndDay)

