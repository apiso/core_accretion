"""
This module calculates the Rosseland mean opacity for the opacity spectra of
D'Alession et al. 2006
"""

import numpy as np
import math
from utils.userpath import userpath
from utils.constants import c, kb, h
from scipy.integrate import odeint
from scipy.interpolate import interp1d

def opacityread(filename):
    lambd = []
    k = []
    f = open(userpath + '/dat/DAlessio_2006_opacities/' + filename, 'r')
    f.readline()
    f.readline()
    for line in f:
        lambd = np.append(lambd, float(line.split()[0]))
        k = np.append(k, float(line.split()[1]))
    f.close()
    return lambd, k


def kR(T, filename, dtg = 0.01):
    lambdm, k = opacityread(filename)
    lambd = lambdm * 10**(-4) #convert from microns to cm
    nu = c / lambd

    dBdT = 2 * h**2 * nu**4 / (c**2 * kb * T**2) * np.exp(h * nu / kb / T) / \
           (np.exp(h * nu / kb / T) - 1)**2

    intdBdTdnu = 8 * kb**4 * np.pi**4 * T**3 / (15 * c**2 * h**3)
        
    kr = 0
    for i in range(len(dBdT) - 1):
        temp = 0.5 * (1 / k[i] * dBdT[i] + 1 / k[i + 1] * dBdT[i + 1]) * \
                   (nu[i] - nu[i + 1])
        if math.isnan(temp):
            temp = 0
        kr = kr + temp
        
    kr = 1./kr * intdBdTdnu
    
    return kr * dtg

def rosseland_opacityread(filename):
    lambd = []
    k = []
    f = open(userpath + '/dat/' + filename, 'r')
    for line in f:
        lambd = np.append(lambd, float(line.split()[0]))
        k = np.append(k, float(line.split()[1]))
    f.close()
    return lambd, k

T, k = rosseland_opacityread('kR_amax_1cm_p_35.txt')

interp_fn = interp1d(T, k)

def interp_opacity(temp):
    return interp_fn(temp)













    
    
