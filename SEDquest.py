from astropy import units as u
from astropy.table import QTable
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt




##############################
# SED generator
# 
# obj = SEDgenerator()
# obj.setcollarea(collarea)
# obj.readcollarea(collarea)
# obj.setsource(sourcespectrum)
# obj.observation(obstime)
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

  def get_collectionarea(self):
    # print("hello")
    return self.CollectionArea
  
  def observation(self,sourcefunction, obstime=60.0*u.min):
    # https://stackoverflow.com/questions/706721/how-do-i-pass-a-method-as-a-parameter-in-python
    # sourcefunction is an object of a function, which is based on a class implemented with a memberfunction 'calc'.
    # minE= 10
    # maxE = 1e5
    # RefEnergies=np.logspace(np.log10(minE),np.log10(maxE))
    
    for e_ref, e_min, e_max, Aeff in zip(self.CollectionArea['e_ref'].value,self.CollectionArea['e_min'].value, self.CollectionArea['e_max'].value, self.CollectionArea['Aeff']):
      flux,_=integrate.quad(sourcefunction, a=e_min, b=e_max)
      # print(flux)
      nevent=flux/u.cm/u.cm/u.s*Aeff*obstime.to(u.s)
      print(nevent,": E=",e_ref,"GeV, flux=",flux,":",Aeff,"->", sourcefunction(e_ref))

    # fluxes=sourcefunction(self.CollectionArea['e_ref'].value)
    # nevents=fluxes*self.CollectionArea['Aeff']*obstime*self.CollectionArea['e_ref'].unit
    # print(self.CollectionArea['e_ref'],nevents)
    # print("****sourcefunction : ",len(nevents))
    # print(self.CollectionArea['e_ref'],self.CollectionArea)
    # print("****collarea ",len(self.CollectionArea))
    

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
    return x*x*1e-6*p0*pow(x/300.,-p1+p2*np.log10(x/300.))
  # TeV/cm2 s

  @staticmethod
  def crabsedmagicnew(x): #Crab MAGIC, arXiv: 1406.6892: Stereo before upgrade
    return x*x*1.e-6*3.23e-11*pow(x/1000.,-2.47-0.24*np.log10(x/1000.))

  @staticmethod
  def crabsedmagicpug(x):#Crab MAGIC, arXiv: 1409.5594: Stereo post upgrade
    return x*x*1.e-6*3.39e-11*pow(x/1000.,-2.51-0.21*np.log10(x/1000.))

  @staticmethod
  def fittedfunc(x,a,b,c,E_0):
  # Function expression: dF/dE = [0]*((x/473.00)^([1]-([2]*([2]*log10(x/473.00)))))
    return x*x*1.e-6*a*pow(x/E_0,b-c*c*np.log10(x/E_0))   
  
#####################################
#  Reference Spectra 
#####################################
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
  def __call__(self, x):
    return self.func(x)
  
class spectrumPL:
  def __init__(self, amplitude,p):
    self.amplitude=amplitude
    self.plindex=p
  def set_amplitude(self, amplitude):
    self.amplitude=amplitude
  def set_plindex(self, plindex):
    self.plindex=plindex
  def __call__(self, x):
    return self.amplitude*pow(x/300., -self.plindex)

class spectrumEPL:
  def __init__(self, amplitude, plindex, Ecut):
    self.amplitude=amplitude
    self.plindex=plindex
    self.cutoffE=Ecut
  def set_plindex(self, plindex):
    self.plindex=plindex
  def set_cutoff(self, Ecut):
    self.cutoffE=Ecut
  def set_amplitude(self, amplitude):
    self.amplitude=amplitude
  def __call__(self, x):
    return self.amplitude*pow(x/300., - self.plindex)

class compositeSpectrum:
  def __init__(self, name):
    self.name = name
    self.functions=[]

  def add_func(self, newfunction):
    self.functions.append(newfunction)

  def __call__(self,x):
    y=np.zeros(len(x))
    for func in self.functions:
      y= y+func(x)
    return y

  

