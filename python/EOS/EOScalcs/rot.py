"""
a module for H2 rotational stat. mech. calculations
>>> import rot
or rename
>>> import rot as r
or import a single fn. or value
>>> from rot import thetaR
for interactive use changes must be reloaded
>>> reload(rot)
"""

import numpy as np
import matplotlib.pyplot as plt

thetaR = 85.5 #rotational temp of H2 in Kelvin
Texample = np.arange(2,300,2.)

def Zrot(T, jmax = 5, verbose = 0):
    """
    partition fn for rotational energy
    """   
    Zr = 0.
    T = np.array(T) #change to numpy array in case list or tuple is passed
    for j in range(jmax + 1):
        weight = (2 - (-1)**j)*(2*j + 1)
        Zrold = Zr
        Zr = Zr + weight*np.exp(-j * (j + 1) * thetaR / T)
        if (verbose) & (j == jmax): print "fractional improvement is ", Zr/Zrold - 1
            
    return Zr

def Urot_on_k(T, jmax = 5, verbose = 0):
    """
    rotational energy (per particle) over k
    """   
    esum = 0.
    T = np.array(T)
    #compute energy sum at a temperature i.e. -(1/kb) * dZ/dbeta
    for j in range(jmax + 1):#print Zr
        weight = (2 - (-1)**j)*(2*j + 1)
        etj = j * (j + 1) * thetaR #rotational energy as a temperature
        esumold = esum
        esum = esum + weight * etj * np.exp(- etj / T)
        if (verbose) & (j == jmax): print "fractional improvement is ", esum/esumold - 1

    return esum / Zrot(T , jmax)

def Cv_on_k(T, jmax = 5, verbose = 0):
    """
    rotational specific heat, per particle over k 
    """   
    esum = 0.
    T = np.array(T)
    #compute squared energy sum i.e. (1/kb**2) * d2Z/dbeta2
    for j in range(jmax + 1):
        weight = (2 - (-1)**j)*(2*j + 1)
        etj = j * (j + 1) * thetaR #rotational energy as a temperature
        esumold = esum
        esum = esum + weight * etj**2 * np.exp(- etj / T)
        if (verbose) & (j == jmax): print "fractional improvement is ", esum/esumold - 1

    #uses Cv/kb = beta**2 * d2lnZ/dbeta2 = beta**2 * ( Z'' / Z - (Z' / Z)**2)
    return 1. / T**2 * (esum / Zrot(T , jmax) - Urot_on_k(T , jmax)**2)

def plotcv(T = Texample, jmax = 5):
    """
    plot rotational cv vs. T
    """
    plt.plot(T, Cv_on_k(T, jmax))
    plt.xlabel('T [K]')
    plt.ylabel('C$_\mathrm{V, rot}$')
    return

