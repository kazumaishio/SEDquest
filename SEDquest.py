from astropy import units as u
from astropy.table import QTable
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt




##############################
# SED generator
# Usage:
# (1) make a class object
# obj = SEDgenerator()
#  
# (2) initial settings 
#   (collection area, background)
# obj.setcollarea(collarea)
# obj.setbackground(bgmodel) 
# obj.readcollarea(collarea)
#  
# (3) generate an observation
# obj.observation(sourcespectrum, obstime)
#  -> proceeds the following steps:
#  - calculate expected # events
#  - generate # events with fluctuation
#  - calculate spectrum with fluctuation  
#  
##############################
# Ingredients: 
# - Spectral shape
#   Flux vs. energy
# e_ref, e_min, e_max, dfde
# - Obs time
#   time
# - Collection area
#  (Migration matrix is omitted)
# e_ref, e_min, e_max, collarea
# 
# Variables:
# - Energy bins (estimated)
#


#####################################
#  SEDgenerator
#####################################
class SEDgenerator:
  def  __init__(self, name):
    self.InstrumentName = name
    self.Observations=QTable()
    self.nobs=0
  # def set_collectionarea(self):
  #   a = np.array([1, 4, 5], dtype=np.int32)
  #   b = [2.0, 5.0, 8.5]
  #   c = ['x', 'y', 'z']
  #   d = [10, 20, 30] * u.m / u.s
  #   self.t = QTable([a, b, c, d],
  #          names=('a', 'b', 'c', 'd'),
  #          meta={'name': 'first table'})
  # def crabsedmagicpug(self,x):#Crab MAGIC, arXiv: 1409.5594: Stereo post upgrade
  #   return x*x*1.e-6*3.39e-11*pow(x/1000.,-2.51-0.21*np.log10(x/1000.))
    
  def set_collectionarea(self,collareapath='./collarea_magic.ecsv'):
    self.CollectionArea= QTable.read(collareapath)
    print('collection area, unit is ',self.CollectionArea['Aeff'].unit)
    self.Observations=self.CollectionArea['e_ref','e_min','e_max' ]
    # print("****collarea ",len(self.CollectionArea))    
    # print(self.CollectionArea['e_ref'],self.CollectionArea)

  def get_collectionarea(self):
    # print("hello")
    return self.CollectionArea
  
  def observation_expected(self,sourcefunction, obstime=60.0*u.min):
    # https://stackoverflow.com/questions/706721/how-do-i-pass-a-method-as-a-parameter-in-python
    # sourcefunction is an object of a function, which is based on a class implemented with a memberfunction 'calc'.
    # minE= 10
    # maxE = 1e5
    # RefEnergies=np.logspace(np.log10(minE),np.log10(maxE))
    # fluxes=sourcefunction(self.CollectionArea['e_ref'].value)
    # nevents=fluxes*self.CollectionArea['Aeff']*obstime*self.CollectionArea['e_ref'].unit
    # print(self.CollectionArea['e_ref'],nevents)
    # print("****sourcefunction : ",len(nevents))
    thisobs=[]
    for e_ref, e_min, e_max, Aeff in zip(self.CollectionArea['e_ref'].value,self.CollectionArea['e_min'].value, self.CollectionArea['e_max'].value, self.CollectionArea['Aeff']):
      if sourcefunction.__class__.__name__=='compositeSpectrum':
        flux=0
        for func in sourcefunction.get_components():
          flux_ind,_=integrate.quad(func, a=e_min, b=e_max)
          flux= flux+flux_ind
      else:
        flux,_=integrate.quad(sourcefunction, a=e_min, b=e_max)
      # print(flux)
      nevent=flux/u.cm/u.cm/u.s*Aeff*obstime.to(u.s)
      (nevent.unit) #for recalculating the unit
      # print(nevent,": E=",e_ref,"GeV, flux=",flux,":",Aeff,"->", sourcefunction(e_ref))
      thisobs.append(nevent)
      
    self.Observations["nevent_expected"]=thisobs
    return thisobs

  def observation(self,sourcefunction, obstime=60.0*u.min):
    # https://stackoverflow.com/questions/706721/how-do-i-pass-a-method-as-a-parameter-in-python
    # sourcefunction is an object of a function, which is based on a class implemented with a memberfunction 'calc'.
    # minE= 10
    # maxE = 1e5
    # RefEnergies=np.logspace(np.log10(minE),np.log10(maxE))
    # fluxes=sourcefunction(self.CollectionArea['e_ref'].value)
    # nevents=fluxes*self.CollectionArea['Aeff']*obstime*self.CollectionArea['e_ref'].unit
    # print(self.CollectionArea['e_ref'],nevents)
    # print("****sourcefunction : ",len(nevents))
    thisobs_nevent=[]
    thisobs_flux=[]
    thisobs_flux_err=[]
    thisobs_sed=[]
    thisobs_sed_err=[]

    for e_ref, e_min, e_max, Aeff in zip(self.CollectionArea['e_ref'].value,self.CollectionArea['e_min'].value, self.CollectionArea['e_max'].value, self.CollectionArea['Aeff']):
      if sourcefunction.__class__.__name__=='compositeSpectrum':
        flux=0
        for func in sourcefunction.get_components():
          flux_ind,_=integrate.quad(func, a=e_min, b=e_max)
          flux= flux+flux_ind
      else :      
        flux,_=integrate.quad(sourcefunction, a=e_min, b=e_max)
      # print(flux)
      nevent=flux/u.cm/u.cm/u.s*Aeff*obstime.to(u.s)
      (nevent.unit) #for recalculating the unit
      nevent_observed=np.random.poisson(nevent)
      # print(nevent_observed,": E=",e_ref,"GeV, flux=",flux,":",Aeff,"->", sourcefunction(e_ref))
      thisobs_nevent.append(nevent_observed)

      flux_observed=nevent_observed/Aeff/obstime.to(u.s)/(e_max-e_min)
      thisobs_flux.append(flux_observed)
      thisobs_sed.append(flux_observed*e_ref*e_ref*1e-6)
      if nevent_observed > 0:
        thisobs_flux_err.append(flux_observed/np.sqrt(nevent_observed))
        thisobs_sed_err.append(flux_observed*e_ref*e_ref*1e-6/np.sqrt(nevent_observed))
      else:
        thisobs_flux_err.append(flux_observed*0)
        thisobs_sed_err.append (flux_observed*e_ref*e_ref*0)

      
    self.Observations["nevent_obs{}".format(self.nobs)]=thisobs_nevent
    self.Observations["flux_obs{}".format(self.nobs)]=thisobs_flux
    self.Observations["flux_err_obs{}".format(self.nobs)]=thisobs_flux_err
    self.Observations["sed_obs{}".format(self.nobs)]=thisobs_sed
    self.Observations["sed_err_obs{}".format(self.nobs)]=thisobs_sed_err
    
    self.nobs=self.nobs+1
    return thisobs_nevent

  def plot_CollectionArea(self):
    plt.plot(self.CollectionArea['e_ref'],self.CollectionArea['Aeff'])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(self.CollectionArea['e_ref'].unit)
    plt.ylabel(self.CollectionArea['Aeff'].unit)  


#####################################
#  Reference Spectra (SED)
#####################################
# minE= 10
# maxE = 1e5
# RefEnergies=np.logspace(np.log10(minE),np.log10(maxE))
# y=LibrefSED.crabsedmagic(RefEnergies)
# print(y)
class LibRefSED:
  @staticmethod
  def crabsedmagic(x): #Crab MAGIC, ApJ 674: Mono
    p0=5.99999999999999996e-10
    p1=2.31000000000000005e+00
    p2=-2.60000000000000009e-01
    return x*x*1e-6*p0*pow(x/300.,-p1+p2*np.log10(x/300.))*u.Unit('TeV/cm2 s')
  # TeV/cm2 s

  @staticmethod
  def crabsedmagicnew(x): #Crab MAGIC, arXiv: 1406.6892: Stereo before upgrade
    return x*x*1.e-6*3.23e-11*pow(x/1000.,-2.47-0.24*np.log10(x/1000.))*u.Unit('TeV/cm2 s')

  @staticmethod
  def crabsedmagicpug(x):#Crab MAGIC, arXiv: 1409.5594: Stereo post upgrade
    return x*x*1.e-6*3.39e-11*pow(x/1000.,-2.51-0.21*np.log10(x/1000.))*u.Unit('TeV/cm2 s')

  @staticmethod
  def fittedfunc(x,a,b,c,E_0):
  # Function expression: dF/dE = [0]*((x/473.00)^([1]-([2]*([2]*log10(x/473.00)))))
    return x*x*1.e-6*a*pow(x/E_0,b-c*c*np.log10(x/E_0))*u.Unit('TeV/cm2 s')
  
#####################################
#  Reference Spectra 
#####################################
# Usage: 
#  register an object of this class to spectrumRef class object. 
#  objspectrum = spectrumRef(LibRefFlux.a_class_method)
class LibRefFlux:
  @staticmethod
  def crabsedmagic(x): #Crab MAGIC, ApJ 674: Mono
    p0=5.99999999999999996e-10
    p1=2.31000000000000005e+00
    p2=-2.60000000000000009e-01
    return p0*pow(x/300.,-p1+p2*np.log10(x/300.))
  # / GeV cm2 s

  @staticmethod
  def crabsedmagicnew(x): #Crab MAGIC, arXiv: 1406.6892: Stereo before upgrade
    return 3.23e-11*pow(x/1000.,-2.47-0.24*np.log10(x/1000.))

  @staticmethod
  def crabsedmagicpug(x):#Crab MAGIC, arXiv: 1409.5594: Stereo post upgrade
    return 3.39e-11*pow(x/1000.,-2.51-0.21*np.log10(x/1000.))

  @staticmethod
  def crabfittedfunc(x,a,b,c,E_0):
  # Function expression: dF/dE = [0]*((x/473.00)^([1]-([2]*([2]*log10(x/473.00)))))
    return a*pow(x/E_0,b-c*c*np.log10(x/E_0))   
  
  @staticmethod
  def mrk501fittedfunc(x):
    p0=4.6194e-10
    p1=-2.20794
    p2=85.7847
    p3=2.31068
    return p0*pow(x/409.00,p1)*np.exp(-pow(x/p2/p2,p3))

  
#####################################
#  Spectral functions
#####################################
class spectrumRef:
  def __init__(self, func):  
    self.func =func
  def calcSED(self,x): #TeV/cm2 s
    return x*x*1e-6*self.func(x)*u.Unit('TeV/cm2 s')
  def __call__(self, x): #/GeV cm2 s
    return self.func(x)
  
class spectrumPL:
  def __init__(self, amplitude,p):
    self.amplitude=amplitude
    self.plindex=p
  def set_amplitude(self, amplitude):
    self.amplitude=amplitude
  def set_plindex(self, plindex):
    self.plindex=plindex
  def calcSED(self,x): #TeV/cm2 s
    return x*x*1e-6*self.amplitude*pow(x/300., -self.plindex)*u.Unit('TeV/cm2 s')
  def __call__(self, x): #/GeV cm2 s
    return self.amplitude*pow(x/300., -self.plindex)

class spectrumEPL:
  def __init__(self, amplitude, plindex, cutoffE,cutCurvature):
    self.amplitude=amplitude
    self.plindex=plindex
    self.cutoffE=cutoffE
    self.cutCurvature=cutCurvature
  def set_plindex(self, plindex):
    self.plindex=plindex
  def set_cutoff(self, cutoffE):
    self.cutoffE=cutoffE
  def set_amplitude(self, amplitude):
    self.amplitude=amplitude
  def calcSED(self,x): #TeV/cm2 s
    return x*x*1e-6*self.amplitude*pow(x/300., - self.plindex)*np.exp(-pow(x/self.cutoffE/self.cutoffE,self.cutCurvature))*u.Unit('TeV/cm2 s')
  def __call__(self, x): #/GeV cm2 s
    return self.amplitude*pow(x/300., - self.plindex)*np.exp(-pow(x/self.cutoffE/self.cutoffE,self.cutCurvature))

class compositeSpectrum:
  def __init__(self, name):
    self.name = name
    self.functions=[]
  def add_func(self, newfunction):
    self.functions.append(newfunction)
  def calcSED(self,x): #TeV/cm2 s
    y=np.zeros(len(x))
    for func in self.functions:
      y= y+x*x*1e-6*func(x)*u.Unit('TeV/cm2 s')
    return y
  def __call__(self,x): #/GeV cm2 s
    # NOTE: this does not work for numpy.integrate.quad function.
    # Thus, each component needs to be processed separately outside this class.
    # -> use get_components() function.
    y=np.zeros(len(x))
    for func in self.functions:
      # print(x)
      # print(func(x))
      y= y+func(x)
    return y
  def get_components(self):
    return self.functions
  

  

