from sys import exit
try:
  import numpy
except ImportError:
  print 'Failure to import numpy module which is critical for this.'
  exit(-1)

class Kernels:
    '''
    Linear kernel models
    '''
    def __init__(self,vza,sza,raa,critical=1,RossHS=False,RecipFlag=True,HB=2.0,BR=1.0,MODISSPARSE=True,MODISDENSE=False,RossType='Thick',normalise=1,normalize=0,LiType='Sparse',doIntegrals=False,BSAangles=[],nbar=0.0):
        '''
        The class creator sets up the kernels for some angle set. Default Li is MODISSPARSE parameter set
	The kernels are accessible from:
		self.Isotropic
		self.Ross
		self.Li
	The angles are accesible from:
		self.vza (or self.vzaDegrees)
		self.sza (or self.szaDegrees)
		self.raa (or self.raaDegrees)
		N.B. Hot spot direction is vza == sza and raa = 0.0
	Kernels integrals are acessible from:
		self.BSAangles (angles in degrees)
		self.BSA_Isotropic (directional-hemispherical integral of self.Isotropic)
		self.BSA_Ross (directional-hemispherical integral of self.Ross)
		self.BSA_Li (directional-hemispherical integral of self.Li)
		self.WSA_Isotropic (bi-hemispherical integral of self.Isotropic)
		self.WSA_Ross (bi-hemispherical integral of self.Ross)
		self.WSA_Li (bi-hemispherical integral of self.Li)
		N.B. You need to set the doIntegrals flag to True on creating an instance of the kernels class if you 
		want access to integrals. The processing takes a bit of time.
	Printing methods are available:
		self.printIntegrals(header=True,reflectance=False)		
		self.printKernels(header=True,reflectance=False)

	Required parameters:

        @param vza: an array containg view zenith angles in degrees
        @param sza: an array containing solar zenith angles in degrees
        @param raa: an array containg relative azimuth angles in degrees

	Options:
        @option critical=1: set to 1 to exit on error, 0 not to
        @option RecipFlag=True: Li reciprocal flag
        @option HB: Li kernel parameter HB 
        @option BR: Li kernel parameter
        @option MODISSPARSE: set to True for default MODIS Li Sparse parameters (overrides BR and HB to 2.0 and 1.0)
        @option MODISDENSE: set to True for default MODIS Li Dense parameters (override BR and HB to 2.0 and 2.5)
        @option RossType: set to 'Thin' for Ross Thin (default) else 'Thick'
        @option LiType: set to 'Sparse' for LiSparse (default). Other options: 'Roujean', 'Dense'
        @option normalise: set to 1 to make kernels 0 at nadir view illumination (default), set to 0 for no normalisation (can also use US spelling, i.e. normalize)
        @option doIntegrals: set to True to calculate integrals of kernels numerically. Set to False not to calculate them. At some point will have Approx flag here as well.
        @option BSAangles: solar zenith angles at which to calculate directional-hemispherical integral of kernels (default 0-89 in steps of 1 degree). Units: degrees.
        @option nbar: the sza at which the isotropic term is set to if normalise=1 is turned on (default 0)

	Notes:
	CRITICAL: requires numpy
 	Requires sys.exit but this is part of any std installation.
	If you do integrals, this also requires scipy (or rather scipy.integrate)
	If you want to plot graphs, requires pylab.
	If you want to mimic the results in Wanner et al. 1995, I've set a special function called self.mimic at the end here. This requires pylab to plot the data.

	Example of use:

	  see mimic()
	
        '''
	self.__setup(critical=critical,RecipFlag=RecipFlag,RossHS=RossHS,HB=HB,BR=BR,MODISSPARSE=MODISSPARSE,MODISDENSE=MODISDENSE,RossType=RossType,normalise=normalise,normalize=normalize,LiType=LiType,doIntegrals=doIntegrals,BSAangles=BSAangles,nbar=nbar)
        self.setAngleInfo(vza,sza,raa)
        self.__doKernels()
        self.__postProcess()
 
    def __setup(self,critical=1,RecipFlag=True,RossHS=False,HB=2.0,BR=1.0,MODISSPARSE=True,MODISDENSE=False,RossType='Thick',normalise=1,normalize=0,LiType='Sparse',doIntegrals=True,BSAangles=[],nbar=0.0):
        self.nbar = nbar
        self.__NEARLYZERO = 1e-20
        self.critical = critical
	self.FILE = -1
	self.outputFileName = 'stdout'
        # kernel options etc.
        self.LiType = LiType
        self.RossHS = RossHS
        self.doIntegrals = doIntegrals
        if MODISDENSE == True:
            LiType = 'Dense'
            self.HB = 2.0
            self.BR = 2.5
        else:
            if MODISSPARSE == True:
                LiType = 'Sparse'
                self.HB = 2.0
                self.BR = 1.0
            else:
                self.HB = HB
                self.BR = BR
        #self.LiType = LiType
        self.RossType = RossType
        self.normalise = normalise
        self.RecipFlag = RecipFlag
        # some useful numbers
        self.__M_PI = numpy.pi
        self.__M_PI_2 = self.__M_PI * 0.5
        self.__M_PI_4 = self.__M_PI * 0.25
        self.__M_1_PI = 1.0/self.__M_PI

        self.normalise = 0
        self.__integrateKernels(BSAangles=BSAangles)

        if (normalise >= 1 or normalize >= 1):
            self.normalise = max(normalise,normalize)

    def __postProcess(self):
	'''
	Private method for dealing with normalisation
	'''
        self.LiNorm = 0.
        self.RossNorm = 0.
        self.IsotropicNorm = 0.
        # if we are normalising the last element of self.Isotropic, self.Ross and self.Li  contain the nadir-nadir kernel
        if self.normalise >= 1:
            # normalise nbar-nadir (so kernel is 0 at nbar-nadir)
            self.RossNorm = self.Ross[-1]
            self.LiNorm = self.Li[-1]
            self.Ross = self.Ross - self.RossNorm
            self.Li  = self.Li - self.LiNorm
            # depreciate length of arrays (well, teh ones we'll use again in any case)
            self.Ross = self.Ross[0:-1]
	    self.Li = self.Li[0:-1]
            self.Isotropic = self.Isotropic[0:-1]
            self.vzaDegrees = self.vzaDegrees[0:-1]
            self.szaDegrees = self.szaDegrees[0:-1]
            self.raaDegrees = self.raaDegrees[0:-1]
            self.N = len(self.vzaDegrees)
            self.vza = self.vza[0:-1]
            self.sza = self.sza[0:-1]
            self.raa = self.raa[0:-1]

    def __doKernels(self):
	'''
	Private method to run the various kernel methods
	'''
        # the kernels
        self.IsotropicKernel()
        self.RossKernel()
        self.LiKernel()

    def setAngleInfo(self,vza,sza,raa):
	'''
	Public method to store and organise the input angle data
	'''
        self.vzaDegrees = numpy.array([vza]).flatten()
        self.szaDegrees = numpy.array([sza]).flatten()
        self.raaDegrees = numpy.array([raa]).flatten()
        self.N = len(self.vzaDegrees)
        
        if(self.N != len(self.szaDegrees) or self.N != len(self.raaDegrees)):
            self.error('kernels: inconsistent number of samples in vza, sza and raa data: ' + str(len(self.vzaDegrees)) + ', ' + str(len(self.szaDegrees)) + ', ' + str(len(self.raaDegrees)),critical=self.critical)
            print self.vzaDegrees
            print self.szaDegrees
            print self.raaDegrees
            return [-1]
        
        if (self.normalise >= 1):
            # calculate nadir term by extending array
            self.vzaDegrees = numpy.array(list(self.vzaDegrees) + [0.0]).flatten()
            self.szaDegrees = numpy.array(list(self.szaDegrees) + [self.nbar]).flatten()
            self.raaDegrees = numpy.array(list(self.raaDegrees) + [0.0]).flatten()
            # not N is one too many now
            self.N = len(self.vzaDegrees)

        self.vza = self.dtor(self.vzaDegrees)
        self.sza = self.dtor(self.szaDegrees) # -1 to make HS direction for raa = 0
        self.raa = self.dtor(self.raaDegrees)
        w = numpy.where(self.vza < 0)[0]
        self.vza[w] = -self.vza[w]
        self.raa[w] = self.raa[w] + numpy.pi
        w = numpy.where(self.sza < 0)[0]
        self.sza[w] = -self.sza[w]
        self.raa[w] = self.raa[w] + numpy.pi


    def __integrateKernels(self,BSAangles=[]):
	'''
	Private method to call integration functions for the kernels


         NB - this overwrites all kernel info ... so be careful how/where you call it
        @option: BSAangles=[] allows the user to set the sza angles at which directional-hemispherical intergal is calculated, else steps of 1 degree from 0 to 89 (though I wouldnt trust it down to 90)
        This function can be rather slow, so using fewer samples or an approximate function may be a god idea
	'''
        if (self.doIntegrals == False):
            return;
        try: 
          import scipy.integrate
          if BSAangles == []:
            BSAangles = numpy.array(range(90))*1.0

          self.BSAangles = numpy.array(BSAangles).flatten()

          # isotropic integral
          self.BSA_Isotropic = numpy.zeros(len(self.BSAangles))+1.0
          self.BSA_Ross = numpy.zeros(len(self.BSAangles))
          self.BSA_Li = numpy.zeros(len(self.BSAangles))
          self.BSA_Isotropic_error = numpy.zeros(len(self.BSAangles))
          self.BSA_Ross_error = numpy.zeros(len(self.BSAangles))
          self.BSA_Li_error = numpy.zeros(len(self.BSAangles))
        
          i = 0
          mu = numpy.cos(self.BSAangles*numpy.pi/180.)
          for sza in self.BSAangles:
            # ross integral
            self.BSA_Ross[i], self.BSA_Ross_error[i] = scipy.integrate.dblquad(RossFunctionForIntegral,0.0, 1.0, __gfun, __hfun, args=(sza,self))
            self.BSA_Li[i], self.BSA_Li_error[i] = scipy.integrate.dblquad(LiFunctionForIntegral,0.0, 1.0, __gfun, __hfun, args=(sza,self))
            i = i + 1
          self.WSA_Ross =  -2.0 * scipy.integrate.simps(self.BSA_Ross * mu,mu)
          self.WSA_Li =  -2.0 * scipy.integrate.simps(self.BSA_Li * mu,mu)
        except ImportError:
          self.error('Warning: failure to import scipy.integrate module: cannot calculate integrals',critical=False)
        return
        
    def __GetPhaang(self):
        '''
        Private method to calculate Phase angle component of kernel
        '''
        self.__cosphaang = self.__cos1*self.__cos2 + self.__sin1*self.__sin2*self.__cos3
        # better check the bounds before arccos ... just to be safe
        w = numpy.where(self.__cosphaang < -1)[0]
        self.__cosphaang[w] = -1.0
        w = numpy.where(self.__cosphaang > 1)[0]
        self.__cosphaang[w] = 1.0       
        self.__phaang = numpy.arccos(self.__cosphaang)
        self.__sinphaang = numpy.sin(self.__phaang)
        return

    def __RossKernelPart(self):
        '''
	Private method to calculate main part of Ross kernel
	'''
        self.__cos1 = numpy.cos(self.vza)
        self.__cos2 = numpy.cos(self.sza)
        
        self.__sin1 = numpy.sin(self.vza)
        self.__sin2 = numpy.sin(self.sza)
        self.__cos3 = numpy.cos(self.raa)
        self.__GetPhaang()
        self.rosselement = (self.__M_PI_2 - self.__phaang)*self.__cosphaang+self.__sinphaang
        return

    def GetDistance(self):
	'''
	Private method to get distance component of Li kernels
	'''
        temp = self.__tan1*self.__tan1+self.__tan2*self.__tan2-2.*self.__tan1*self.__tan2*self.__cos3;
        w = numpy.where(temp < 0)[0]
        temp[w] = 0.0
        self.__temp = temp # used by other functions ??
        distance = numpy.sqrt(temp)
        return distance

    def GetpAngles(self, tan1):
        '''
        Private method to do B/R transformation for ellipse shape
        '''
        t = self.BR * tan1
        w = numpy.where( t < 0.)[0]
        t[w] = 0.0
        angp = numpy.arctan(t)
        s = numpy.sin(angp)
        c = numpy.cos(angp)
        # have to make sure c isnt 0
        w = numpy.where(c == 0)[0]
        c[w] = self.__NEARLYZERO
        return c,s,t

    def GetOverlap(self):
        '''
        Private method to do HB ratio transformation
        '''
        self.__temp = 1./self.__cos1 + 1./self.__cos2

        self.__cost =  self.HB * numpy.sqrt(self.__distance * self.__distance + self.__tan1 * self.__tan1 * self.__tan2 * self.__tan2 * self.__sin3 * self.__sin3) / self.__temp;
        w = numpy.where(self.__cost < -1)[0]
        self.__cost[w] = -1.0
        w = numpy.where(self.__cost > 1.0)[0]
        self.__cost[w] = 1.0
        self.__tvar = numpy.arccos(self.__cost)
        self.__sint = numpy.sin(self.__tvar)
        self.__overlap = self.__M_1_PI * (self.__tvar - self.__sint * self.__cost) * self.__temp
        w = numpy.where(self.__overlap < 0)[0]
        self.__overlap[w] = 0.0
        return
         
    def RoujeanKernel(self):
        '''
        Private method - call to calculate Roujean shadowing kernel
        '''
        # first make sure its in range 0 to 2 pi
        self.__phi = numpy.abs((self.raa % (2.*numpy.pi)))
        self.__cos3 = numpy.cos(self.__phi)
        self.__sin3 = numpy.sin(self.__phi)
        self.__tan1 = numpy.tan(self.sza)
        self.__tan2 = numpy.tan(self.vza)

        self.__distance = self.GetDistance()
        self.Li = 0.5 * self.__M_1_PI * ((self.__M_PI - self.__phi) * self.__cos3 + self.__sin3) * self.__tan1 * self.__tan2 - self.__M_1_PI * (self.__tan1 + self.__tan2 + self.__distance);
        return

    def LiKernel(self):
        '''
        Private method - call to calculate Li Kernel
        '''
        # at some point add in LiGround kernel & LiTransit
        if self.LiType == 'Roujean':
            return self.RoujeanKernel()
        # first make sure its in range 0 to 2 pi
        self.__phi = numpy.abs((self.raa % (2.*numpy.pi)))
        self.__cos3 = numpy.cos(self.__phi)
        self.__sin3 = numpy.sin(self.__phi)
        self.__tanti = numpy.tan(self.sza)
        self.__tantv = numpy.tan(self.vza)
        self.__cos1, self.__sin1, self.__tan1 = self.GetpAngles(self.__tantv);
        self.__cos2, self.__sin2, self.__tan2 = self.GetpAngles(self.__tanti);
        self.__GetPhaang(); # sets cos & sin phase angle terms 
	self.__distance = self.GetDistance(); # sets self.temp
	self.GetOverlap(); # also sets self.temp
        if self.LiType == 'Sparse':
            if self.RecipFlag == True:
                self.Li = self.__overlap - self.__temp + 0.5 * (1. + self.__cosphaang) / self.__cos1 / self.__cos2;
            else:
                self.Li = self.__overlap - self.__temp + 0.5 * (1. + self.__cosphaang) / self.__cos1;
        else:
 	    if self.LiType == 'Dense':
                if self.RecipFlag:
                    self.Li = (1.0 + self.__cosphaang) / (self.__cos1 * self.__cos2 * (self.__temp - self.__overlap)) - 2.0;
                else:
                    self.Li = (1.0 + self.__cosphaang) / (self.__cos1 * (self.__temp - self.__overlap)) - 2.0;
            else:
	        B = self.__temp - self.__overlap
	        w = numpy.where(B <= 2.0)
	        self.Li = B*0.0
                if self.RecipFlag == True:
	            Li = self.__overlap - self.__temp + 0.5 * (1. + self.__cosphaang) / self.__cos1 / self.__cos2;
	        else:
		    Li = self.__overlap - self.__temp + 0.5 * (1. + self.__cosphaang) / self.__cos1;
	        self.Li[w] = Li[w]

		w = numpy.where(B > 2.0)
		if self.RecipFlag:
		    Li = (1.0 + self.__cosphaang) / (self.__cos1 * self.__cos2 * (self.__temp - self.__overlap)) - 2.0;
		else:
		    Li = (1.0 + self.__cosphaang) / (self.__cos1 * (self.__temp - self.__overlap)) - 2.0;
		self.Li[w] = Li[w]
        return
        
    def IsotropicKernel(self):
        '''
        Public method - call to calculate Isotropic kernel
        '''
        # default behaviour
        self.Isotropic = numpy.zeros(self.N)+1.0
        return

    def RossThin(self):
	'''
	Public method - call to calculate RossThin kernel
	'''
        self.__RossKernelPart()
        self.rosselement = self.rosselement/(self.__cos1*self.__cos2)
        return;

    def RossThick(self):
	'''
	Public method - call to calculate RossThick kernel
	'''
        self.__RossKernelPart()
        self.rosselement = self.rosselement/(self.__cos1+self.__cos2)
        return;

    def RossKernel(self):
        '''
        Public method - call to calculate Ross Kernel

	Sets self.Ross to result
        '''
        if self.RossType == 'Thin':
            self.RossThin()
        else:
            self.RossThick()
        self.Ross = self.rosselement
	if self.RossHS != False:
	    if self.RossHS == True:
		self.RossHS = 0.25
            self.Ross = self.Ross * (1 + 1/(1 + self.__phaang/self.RossHS))
    

    def dtor(self,x):
	'''
	Public method to convert degrees to radians

        arguments:
          x single number or array of numbers in degrees

        returns:
          float point cast of x, converted to radians

	'''
        return x*numpy.pi/180.0

    def rtod(self,x):
	'''
	Public method to convert radians to degrees

	arguments:
	  x single number or array of numbers in radians

	returns:
	  float point cast of x, converted to degrees
	'''
        return x*180./numpy.pi
    
    def error(self,msg,critical=0,newline=1,code=-1):
        '''
        Public method to do Class error reporting
        @param msg: error message
        @param critical: set to 1 if require exit (default critical=0)
        @param newline: set to 0 if newline not required (default newline=0)
        @param code: error code reported on exit if critical error (default code=-1)

	returns:
	  error code (default -1)

        '''
        if newline == 1:
            nl = '\n'
        else:
            nl = ''
        print msg + nl
        if critical == 1:
          exit([code])
        return code

    def printIntegrals(self,header=True,reflectance=False):
	'''
	Public method to print kernel integrals

	By default this will go to stdout
	To output to some file Name, call self.outputFile(file=Name) prior to calling this method

        returns:
	  nothing

	'''
        if(header == True):
            self.printer('# ' + str(self.N) + ' samples Ross: ' + self.RossType + ' Li: ' + self.LiType + ' Reciprocal: ' + str(self.RecipFlag) + ' normalisation: ' + str(self.normalise) + ' HB ' + str(self.HB) + ' BR ' + str(self.BR) + '\n');
            self.printer('# WSA: Isotropic 1.0 Ross ' + str(self.WSA_Ross) + ' Li ' + str(self.WSA_Li))
            self.printer('# 1: SZA (degrees) 2: BSA Isotropic 3: BSA Ross 4: BSA Li')
            if (reflectance == True):
                self.printer(' ');
            self.printer('\n');
            
        for i in range(len(self.BSAangles)):
            self.printer(str(self.BSAangles[i]) +  ' ' + str(self.BSA_Isotropic[i]) + ' ' + str(self.BSA_Ross[i]) + ' ' + str(self.BSA_Li[i]))
            # print refl data if wanted
            if (reflectance == True):
                self.printer(' ');
            self.printer('\n');
        return

    def printKernels(self,header=True,reflectance=False,file=False):
        '''
        Public utility method to print kernel values  

	The method has no arguments  

	It has the following options:
	   header=True|False         (default=True) : write two lines of metedata at the head of the file
	   reflectance=True|False    (default=False): pass through some reflectance data to prtint into the file (not implemented)
	   file=False|Name           (default=False): If set to false write to write to whatever file is open (self.outputFile())
			                              If Name is set then write data to this filename and set self.outputFile(file=Name)
	returns:
	  nothing

	'''
        if(file != False):
	    if(file != self.outputFileName):
		self.outputFile(file=file)

        if(header == True):
            self.printer('# ' + str(self.N) + ' samples Ross: ' + self.RossType + ' Li: ' + self.LiType + ' Reciprocal: ' + str(self.RecipFlag) + ' normalisation: ' + str(self.normalise) + ' HB ' + str(self.HB) + ' BR ' + str(self.BR) + '\n');
            self.printer('# 1: VZA (degrees) 2: SZA (degrees) 3: RAA (degrees) 4: Isotropic 5: Ross 6: Li')
            if (reflectance == True):
                self.printer(' ');
            self.printer('\n');
            
        for i in range(self.N):
            self.printer(str(self.vzaDegrees[i]) + ' ' + str(self.szaDegrees[i]) + ' ' + str(self.raaDegrees[i]) + ' ' + str(self.Isotropic[i]) + ' ' + str(self.Ross[i]) + ' ' + str(self.Li[i]))
            # print refl data if wanted
            if (reflectance == True):
                self.printer(' ');
            self.printer('\n');
        return

    def outputFile(self,file=False,close=True,mode="w"):
	'''
	Public methods to return or set output file name

	No arguments

	Options:
	file=False|Name (default False) : if file not False, set the output file name to Name and open it (closing the previous output file unless close=False is set)
	close=True|False (default True) : close any previous file before opening this one
	mode="w"|"a"    (default "w")   : specify whether opening file for creation or appending

        returns:
          string of filename or False if file could not be opened
	'''
        if (file != False):
	  if (close == True and self.FILE >= 0 and self.FILE != False):
            try:
              self.FILE.close()
            except IOError:
              self.error('failure to close file ' + self.outputFileName)
	  self.outputFileName=file
	  if (file == "stdout"):
	    self.FILE = -1
	  elif (file == "stderr"):
	    self.FILE = -2
	  else:
	    try:
	      self.FILE = open(self.outputFileName,mode)
            except IOError:
	     self.error('failure to open file '  + self.outputFileName + ' for writing, mode ' + mode)
	     self.FILE = False
	     return False
	     
	return self.outputFileName

    def printer(self,msg):
        '''
        Public print method 

	arguments: msg (string)

        write ascii string msg to stdout (default) or file self.outputFile())
        returns:
          nothing
        ''' 
        if (self.FILE == -1):
	  print msg,
	else:
	  self.FILE.write(msg)


# some things required for the numerical integration

def _Kernels__gfun(x):
    return 0.0

def _Kernels__hfun(x):
    return 2.0*numpy.pi

def RossFunctionForIntegral(phi,mu,sza,self):
    #print phi
    #print mu
    #print sza
    #print '========'
    vza = numpy.arccos(mu)
    raa = self.rtod(phi)
    self.setAngleInfo(vza,sza,raa)
    self.RossKernel()
    return mu * self.Ross[0] / numpy.pi

def LiFunctionForIntegral(phi,mu,sza,self):
    #print phi
    #print mu
    #print sza
    #print '========'
    vza = numpy.arccos(mu)
    raa = self.rtod(phi)
    self.setAngleInfo(vza,sza,raa)
    self.LiKernel()
    return mu * self.Li[0] / numpy.pi

# test function
def mimic(doPrint=False,doPlot=False,RossHS=False,RecipFlag=False):
    '''
    A test method to reproduce the results in Wanner et al. 1995.

    It provides access to data and plots of various kernels that should match those
    in that publication (see below for flags) or that can be used to visualise other forms of the kernels.

    There are no parameters but the following options:
            doPrint=True|False    : print results to stdout (default doPrint=False)
	    doPlot=True|False     : plot the data in pylab
	    RossHS=True|False     : set the Ross Host Spot term (default=False)
	    RecipFlag=True|False  : set the Li kernels to their reciprocal versions (default=False)
            NoReturn=True|False   : do|don't return anything from the fn (default=False)
  
    Note that to match the graphs in Wanner et al. 1995 you should set RossHS=False,RecipFlag=False
    Note also that for the kernels used in the MODIS algorithm you should set RossHS=False,RecipFlag=True

    The method returns (unless NoReturn==True):
	VZA,SZA,RAA,RossThick,RossThin,LiSparse,LiDense,Roujean,LiTransit
    where all are numy arrays of dimensions 3 x nSamples 
    so:
        VZA[0,:],RossThick[0,:] are the results for sza = 0.0
        VZA[1,:],RossThick[1,:] are the results for sza = 30.0
        VZA[2,:],RossThick[2,:] are the results for sza = 60.0

    and/or it provides data in files/plots of these kernels.
    Note that LiTransit is not included in Wanner et al. 1995.
    '''
    # set up the angles
    r = 89 # do results for +/- r degrees)
    SZAS = numpy.array([0.0,-30.0,-60.0])  # sza
    vza = numpy.array(range(2*r+1))*1.0 - r
    # set up storage info
    RossThick = numpy.zeros([3,len(vza)])
    RossThin = numpy.zeros([3,len(vza)])
    LiSparse = numpy.zeros([3,len(vza)])
    LiDense = numpy.zeros([3,len(vza)])
    Roujean = numpy.zeros([3,len(vza)])
    LiTransit = numpy.zeros([3,len(vza)])
    SZA = numpy.zeros([3,len(vza)])
    VZA = numpy.zeros([3,len(vza)])
    RAA = numpy.zeros([3,len(vza)])
    # fill the angle info
    RossHS=RossHS
    for i in range(len(SZAS)):
  	SZA[i,:] = SZAS[i]
	VZA[i,:] = vza[:] 
	RAA[i,:] = 0.0
        # do the kernels
        kk = Kernels(VZA[i,:] ,SZA[i,:],RAA[i,:],RossHS=RossHS,MODISSPARSE=True,RecipFlag=RecipFlag,normalise=1,doIntegrals=False,LiType='Dense',RossType='Thick')
    	RossThick[i,:] = kk.Ross[:]
        LiDense[i,:] = kk.Li[:]
 	if doPrint == True:
	  kk.printKernels(file='RossThickLiDense.' + str(SZAS[i]) + '.dat')
	  kk.printer('')
	kk = Kernels(VZA[i,:] ,SZA[i,:],RAA[i,:],RossHS=RossHS,MODISSPARSE=True,RecipFlag=RecipFlag,normalise=1,doIntegrals=False,LiType='Sparse',RossType='Thin')
	RossThin[i,:] = kk.Ross[:]
        LiSparse[i,:] = kk.Li[:]
        if doPrint == True:
                kk.printKernels(file='RossThinLiSparse.' + str(SZAS[i]) + '.dat')
                kk.printer('')
	kk = Kernels(VZA[i,:] ,SZA[i,:],RAA[i,:],RossHS=RossHS,MODISSPARSE=True,RecipFlag=RecipFlag,normalise=1,doIntegrals=False,LiType='Roujean',RossType='Thin')
	Roujean[i,:] = kk.Li[:]
	if doPrint == True:
                kk.printKernels(file='RossThinRoujean.' + str(SZAS[i]) + '.dat')
                kk.printer('')
        kk = Kernels(VZA[i,:] ,SZA[i,:],RAA[i,:],RossHS=RossHS,MODISSPARSE=True,RecipFlag=RecipFlag,normalise=1,doIntegrals=False,LiType='Transit',RossType='Thin')
        LiTransit[i,:] = kk.Li[:]
        if doPrint == True:
                kk.printKernels(file='RossThinLiTransit.' + str(SZAS[i]) + '.dat')
                kk.printer('')
    if (doPlot == True):
	try: 
          import pylab
	  x = [-90.0,90.0]
	  y = [0.0,0.0]
	  for i in range(len(SZAS)):
	    sza = SZAS[i]
 	    pylab.clf()
	    pylab.xlabel('View Zenith Angle')
	    pylab.ylabel('Kernel Value')
	    pylab.title('Solar Zenith Angle ' + str(sza) + ' Degrees')
	    pylab.plot(x,y)
	    pylab.plot(kk.vzaDegrees,RossThick[i,:],label='RThick')
            #pylab.plot(kk.vzaDegrees,RossThin[i,:],label='RThin')
            pylab.plot(kk.vzaDegrees,LiSparse[i,:],label='LiSp')
            #pylab.plot(kk.vzaDegrees,LiDense[i,:],label='LiDen')
            #pylab.plot(kk.vzaDegrees,Roujean[i,:],label='Roujean')
            #pylab.plot(kk.vzaDegrees,LiTransit[i,:],label='LiTrans')
	    pylab.axis([-90.0,90.0,-3.0,3.0])
            pylab.legend(loc=0)
            pylab.savefig(str(sza) + '.png')
	    #pylab.show()
        except ImportError:
          kk.error('Error attempting to import pylab ... you cant do any plotting',critical=False)


    return VZA,SZA,RAA,RossThick,RossThin,LiSparse,LiDense,Roujean,LiTransit
