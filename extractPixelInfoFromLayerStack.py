#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

"""
Extract temporal profile data for a single pixel
Calculate local min and max
Compute change points

Authors:   Gerardo Lopez-Saldana <GerardoLopez@isa.utl.pt>
"""
import os

try:
  import numpy
  import matplotlib as mpl
  mpl.use('Agg')
  from pylab import *
except ImportError:
  print 'Numpy/MatPlotLib is not installed'
  exit(-1)

try:
  import osgeo.gdal as gdal
  from osgeo.gdalconst import *
  gdal.UseExceptions()
except ImportError:
  print 'GDAL is not installed.'
  exit(-1)

try:
  sys.path.append('/Users/glopez/GCII/src/ChangePoint')
  import ChangePoint
except ImportError:
  print 'Change point detection source is not in the right path.'
  exit(-1)


def Interpolate(X, Y):
    from scipy.interpolate import interp1d
    # Interpolate where there is no data
    indices = numpy.where(Y<>0.0)
    x_with_data = X[indices[0]]
    y = Y[indices[0]]

    # find the first and last values on X
    index_first_element = numpy.where(X == x_with_data[0])
    index_last_element = numpy.where(X == x_with_data[-1])
    new_x = X[index_first_element[0][0]:index_last_element[0][0]]

    f = interp1d(x_with_data, y)
    interpolated_data = f(new_x)

    return new_x, interpolated_data

def GetPixel(File, xmin=1, ymin=1, xmax=1, ymax=1):
    import sys
    import glob

    try:
        filename = glob.glob(File)
        dataset = gdal.Open( filename[0], GA_ReadOnly )
        BandCount = dataset.RasterCount
    except:
        print "Error:", sys.exc_info()[0]
        exit(-1)

    #print "Opening file", filename[0]

    # Load data
    data = numpy.zeros(BandCount, numpy.float32)
    DoY = numpy.zeros(BandCount, numpy.float32)

    #from IPython import embed
    #ipshell = embed()

    j = 1
    # The additionall loop is due to a some memory problem with GDAL, therefore extract only 30 bands and reset the dataset
    for i in range(BandCount):
        ##try:
        if j <= 30:
            data[i] = dataset.GetRasterBand(i+1).ReadAsArray(xmin-1, ymin-1, 1, 1)
            j += 1
        else:
            dataset = None
            dataset = gdal.Open( filename[0], GA_ReadOnly )
            data[i] = dataset.GetRasterBand(i+1).ReadAsArray(xmin-1, ymin-1, 1, 1)
            j = 1
        ##except RuntimeError:
        ##    print "Error:", sys.exc_info()[0]
        ##    dataset = None
        ##    dataset = gdal.Open( filename[0], GA_ReadOnly )
        ##    data[i] = dataset.GetRasterBand(i+1).ReadAsArray(xmin-1, ymin-1, 1, 1)

    y = data

    dataset = None
    return y

import sys

# Layerstack of BRDF adjusted BB NIR BHR 
File = '/Users/glopez/GCII/data/LTDR/Filtered/1998/NIR/BHR_NIR.1998.vrt'

sample = int(sys.argv[1])
line = int(sys.argv[2])

NumberOfFires = len(sys.argv[3].split(','))
strFires = sys.argv[3].split(',')
FireDoY = numpy.zeros(NumberOfFires)
for i in range(NumberOfFires):
	FireDoY[i] = int(strFires[i])

NumberOfFiresInCluster = sys.argv[4]

xmin = xmax = sample
ymin = ymax = line

y = GetPixel(File, xmin, ymin, xmax, ymax)
x = numpy.array(range(1,365,8))

ax = subplot(111)

# Draw a vertical line at the inital DoY of each month
for month in [1,32,60,91,121,152,182,213,244,274,305,335]:
	l = axvline(x=month, color='black', ls='-.')

# Plot a vertical line indicating the fire DoY
#l = axvline(x=FireDoY, color='red', ls='--')
for fire in FireDoY:
	l = axvline(x=fire, color='red', lw=1.3, label='Active fire')

plot(x, y, 'bx')
plot(x, y, 'b', label="BB BRDF-corrected", lw=0.7)

# Daily broadband NIR non-BRDF-corrected reflectances
DailyReflectances = '/Users/glopez/GCII/data/LTDR/Filtered/1998/DailyNIR/Daily_NIR.1998.vrt'
y_daily = GetPixel(DailyReflectances, xmin, ymin, xmax, ymax)
x_daily = numpy.array(range(1,366))

indices = numpy.where(y_daily==0)
plot(x_daily[indices], y_daily[indices], 'o', color='#C0C0C0')

indices = numpy.where(y_daily<>0)
plot(x_daily[indices], y_daily[indices], 'x', color='#C0C0C0')
plot(x_daily[indices], y_daily[indices], color='#C0C0C0', lw=0.8, label="Daily BB non-BRDF-corrected")

# Daily narrowband NIR non-BRDF-corrected reflectances
##DailyNarrowBandReflectances = '/Users/glopez/GCII/data/LTDR/N14/1998/AVH09C1.NIR.1998.vrt'
##y_daily_narrowband = GetPixel(DailyNarrowBandReflectances, xmin, ymin, xmax, ymax)
##y_daily_narrowband = y_daily_narrowband * 0.0001

##indices = numpy.where(y_daily_narrowband<>0)
##plot(x_daily[indices], y_daily_narrowband[indices], 'x', color='#666666')
##plot(x_daily[indices], y_daily_narrowband[indices], color='#666666', label="Daily NB non-BRDF-corrected")

# The "+1" is mandatory since "diff" reduces the original index number.
x_interp, y_interp = Interpolate(x_daily, y_daily)
indices = numpy.where(y_interp<>0)

#from IPython import embed
#ipshell = embed()

LocalMin = (numpy.diff(numpy.sign(numpy.diff(y_interp[indices]))) > 0).nonzero()[0] + 1 
# Add the last point to the LocalMin time series
#LocalMin = numpy.append(LocalMin, indices[-1])
LocalMax = (numpy.diff(numpy.sign(numpy.diff(y_interp[indices]))) < 0).nonzero()[0] + 1
plot(x_interp[LocalMin], y_interp[LocalMin], 'black', lw=0.7,)
plot(x_interp[LocalMax], y_interp[LocalMax], 'g', lw=0.5,)

# Default call: ChangePoint(data, confidence=95., iterations=1000)
#points = sorted(set(ChangePoint.ChangePoint(y_interp[LocalMin], confidence=99.)))
#points = ChangePoint.ChangePoint(y_interp[LocalMin], confidence=99.)
points = ChangePoint.ChangePoint(y_interp, confidence=80.)
#points = ChangePoint.ChangePoint(y_interp[LocalMin], confidence=80.)
print points
print
#for i, a in zip(points, numpy.split(y_interp, points)):
#    if len(a) == 0:
#        continue
#    print x[i], "==>", a

#for cp in points:
if points:
    cp = points[0]
    confidence = points[1]
    print cp
    index = numpy.where( y_interp == y_interp[cp]  )
    #index = numpy.where( y_interp == y_interp[LocalMin][cp]  )
    l = axvline(x=x_daily[index[0]], lw=0.8, color='#FFA500', label='Change point DoY ' + str(int(x_daily[index[0]])) + ' ' + 'Confidence ' + str(confidence))

# Change font size
import matplotlib.font_manager as fm
properties = fm.FontProperties(size=10)
ax.legend(prop=properties)
#File = sys.argv[1]
axis([1.0,365,0.0,y.max()+0.1])

# Save the plot
filename = 'Sample_' + str(sample) + '_Line_' + str(line) + '_NumberOfFiresInCluster_' + NumberOfFiresInCluster + '.png'
savefig(filename, dpi=140)
close()


