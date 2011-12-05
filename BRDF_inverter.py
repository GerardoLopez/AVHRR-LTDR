# -*- coding: utf-8 -*-
"""
A class to invert surface reflectance using the Lewis & Gomez-Dans algorithm.

:author: Jose Gomez-Dans <j.gomez-dans@geog.ucl.ac.uk>

"""
#import pdb
import numpy
from scipy.linalg import lstsq, inv

from Kernels import Kernels

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class NotEnoughSamplesError(Error):
    """
       Exception raised when not enough samples are available.
       
    Attributes:
    -----------
    :samples: Number of available samples for inversion
    :author: Jose Gomez-Dans <j.gomez-dans@geog.ucl.ac.uk>
    """

    def __init__(self, samples ):
        self.samples = samples

class AngularNormalisation:
    def __init__ ( self, year, win_len=16, \
                  angular_inversion="TemporalWindow", \
                  doIntegrals=False, RossHS=False, RossType='Thick', \
                  LiType='Sparse', normalise=1, \
                  RecipFlag=True, MODISSPARSE=True, sensor="MODIS"  ):

       """
          The class constructor. It sets up loads of things, some of
          which may be overridden
          by other member methods.
        """
        # Define some useful global variables
        self.angular_inversion = angular_inversion
        #this is usueful for DoY translations and so on
        self.year = year
        """The type of angular inverstion required: either temporal
        window or use of cubic temporal function."""
        self.doIntegrals = doIntegrals
        self.RossHS = RossHS
        self.RossType = RossType
        self.LiType = LiType
        self.normalise = normalise
        self.RecipFlag = RecipFlag
        self.MODISSPARSE = MODISSPARSE
        self.win_len = win_len
        # Initially, only MODIS considered, but could add MERIS or others
        # Problem is whether the band uncertainty is known for other
        # sensors.
        self.sensor=sensor
        """The sensor type. Only MODIS considered for now!"""
        if self.sensor == "MODIS":
            self.bu = numpy.array([0.004, 0.015, 0.003, 0.004, 0.013, \
                        0.010, 0.006])
            """Band uncertainties."""
            self.QA_OK = numpy.array([8, 72, 136, 200, 1032, 1288, 2056, \
                        2120, 2184, 2248])
            """Acceptable QA flags."""
            self.wavelengths = numpy.array([645., 858.5, 469., 555., \
                        1240., 1640., 2130.])
            self.nBands = self.wavelengths.shape[0]

    def do_normalisation ( qa, vza, sza, raa, refl, doy, dos,\
                              pasar=None ):
        """A method to carry out angular normalisation of surface
        reflectance data.

        The method checks for a minimum value of available observations (once
        QA has been evaluated), and will rise an exception
        :exc:`NotEnoughSamplesError`

        :param qa: QA flags. Usually, from MODIS, you get int16 values.
        The ones that are acceptable are given in :attr:`self.QA_OK`.
        :param vza: View Zenith Angle in degrees \*100 (as in the
        MOD09/MYD09 products)
        :param sza: Sun Zenith Angle in degrees \*100 (as in the
        MOD09/MYD09 products)
        :param raa: Relative Azimuth Angle in degrees \*100 (as in the
        MOD09/MYD09 products)
        :param refl: The reflectance (\*10000), an array of
        :math:`(7 \times N_{samples})`
        :param doy: An array of days of year in MODIS format (YearDoY).
        We can use fractions to indicate several acquisitions on the same
        date.
        :param dos: The "day of estimation Same format as doy
        :param pasar: If the QA data is already a 0/1 mask, you can pass
        the 0/1 mask as the QA here. Remember to pass it as :obj:`qa` too.
        :returns: :math:`\\rho_{s}`

        """
        if pasar is None:
            pasar = numpy.logical_or.reduce([qa==x for x in self.QA_OK])
        pasar = (pasar==1)
        if(pasar.sum()<=18):  #Not enough samples
            raise NotEnoughSamplesError( pasar.sum() )
            #return None
        vza = vza[pasar]/100.
        sza = sza[pasar]/100.
        raa = raa[pasar]/100.
        refl = refl[:, pasar]/10000.
        doy = doy[pasar]
        doy = doy-self.year*1000

        kk = Kernels(vza, sza, raa, \
                    doIntegrals=self.doIntegrals, RossHS=self.RossHS, \
                    RossType=self.RossType, LiType=self.LiType,\
                    normalise=self.normalise, RecipFlag=self.RecipFlag, \
                    MODISSPARSE=self.MODISSPARSE )

        
        calc_kernels, residuals, uncertainty = self.InvertTemporalWindow \
                    ( kk, refl, doy, dos )

        return calc_kernels, residuals, uncertainty

    def InvertTemporalWindow ( self, kk, refl, doy, dos):
        """Temporal window inversion (AKA "16 day window inversion")

           This method implements the MODIS MOD43 16-day window inversion
           for reflectance values. First, it checks
           that we have enough samples, and if not, raises
           :exc:`NotEnoughSamplesError`.

           :param kk: An instance of the ``Kernel`` object.
           :param refl: The measured reflectance.
           :param doy: Corresponding dates to the reflectance & kernel
           observations.
           :param dos: Day of burn.
           :param win_len: Lenght of the temporal window to invert. Set by
           default to 16 days.
           :returns: :math:`\\rho_{s}`
        """
        #Calculate valid samples to use.
        passer = numpy.logical_and ( doy<=(dos), \
                    doy>(dos-(self.win_len+1)))

        if (passer.sum()<7) :
            #Not enough samples
            raise NotEnoughSamplesError( passer.sum() )
        
        K = numpy.ones([3, passer.sum()])
        K[1, :] = kk.Ross[passer]
        K[2, :] = kk.Li[passer]
        (P, rho_residuals, rank, svals) = lstsq ( K.T, \
                                              refl[:, passer].T)
        # Uncertainty for the f0 paramter is given by
        # :ref:`Lucht & Lewis`
        # s^2 =  MSE * M-1[0,0], where M^{-1} is
        # K*K^{T}
        M = numpy.dot ( K, K.T )
        unc_factor = numpy.linalg.inv ( M ) [0,0]
        unc = rho_residuals * unc_factor
        return P, rho_residuals, unc
