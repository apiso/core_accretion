"""

Constants for atmosphere modelling ; cgs unless otherwise noted

"""
import numpy as np
from collections import namedtuple

#physical
G = 6.673e-08
c = 2.99792458e10
h = 6.626068e-27
kb = 1.3807e-16
mp = 1.67262171e-24
mH = mp
NA = 6.02214129e23
Rb = 83144720.0 #NA*kb #ideal gas constant (per mol not per H!)
RH = Rb/2
RHe = Rb/4
sigma =  5.670373e-5

#unit conversion
cmperau = 14959788700000.0
AU = 1.49597887e13
year = 3.16e7
Myr = 3.16e13
barP = 1e6 #dyne/cm*2 to bar (not that `bar` in matplotlib is problem in pylab)

#astronomical
Me = 5.9742e27
Re = 6.371e8
Msun = 1.9891e33

def RHill(mp, a, Ms = Msun):
    """
    Hill radius in cgs with SMA (`a`) in AU and planet mass (`mp`) also
    in cgs
    """
    return (mp / 3 / Ms)**(1./3) * a * AU

#MMSN disk model
#(could/ should make the deviations from standard MMSN optional or keyword parameters)
def Tdisk(a, FT):  
    """Disk temperature as function of semi-major axis (a in AU) and Ft"""
    return 120 * FT * (a)**(-3./7)

def Pdisk(a, mstar, FSigma, FT): 
    """
    Disk pressure (in dyn cm^-2) as function of semi-major axis (a in AU),
    stellar mass (mstar in solar masses), and factors Fsigma and Ft
    """
    return 1.1 * 10**(-4) * FSigma * np.sqrt(FT * mstar) * a**(-45./14) * 10**6

def Hdisk(a, FT):  #a in AU
    """Disk scaleheight (in cm) as function of semi-major axis and Ft"""
    return 0.021784 * AU * np.sqrt(FT) * (a)**(9./7) 

#opacity laws

def kdust(T, P = None, b = 2, fsolar = 1.):
    """Opacity (generally applicable for beta = 2)"""
    return 2 * fsolar * (T/100.)**b

def kdust10(T, P = None, b = 2, fsolar = 1.):
    """Opacity reduced by a factor of 10"""
    return 2 * fsolar * (T/100.)**b / 10

def kdust100(T, P = None, b = 2, fsolar = 1.):
    """Opacity reduced by a factor of 100"""
    return 2 * fsolar * (T/100.)**b / 100

def kfixed(T = None, P = None, const = 0.01):
    """a silly function but useful if kappa functions expects a temperature."""
    return const

def gammafn(delad):
    """Polytropic index gamma as function of delad."""
    return 1. /(1 - delad) 

#R, Cv, molecular weight for a given H/He mixture fraction Y

def mufn(Y):
    return 4./(2-Y)

def Rfn(Y):
    return Rb / mufn(Y)

def Cvfn(Y, delad):
    """
    Assumes that the He fraction only affects mu, but not delad.
    """
    return Rfn(Y) / (gammafn(delad) - 1)

#tuple of gas, disk and core parameters
params = namedtuple( \
    'params', 'Mco, rco, a, delad, Y, gamma, R, Cv, Pd, Td, kappa')




