"""
This module contains functions that make generating plots easier. Used by
scripts & ipython notebooks in python/Atmscrips.
"""


import numpy
from utils.constants import Me, Re, gammafn, Rfn, Cvfn, Pdisk, Tdisk, kdust, \
     params
from utils.parameters import FSigma, FT, mstar
import matplotlib.pyplot as plt

from RadSGPoly.profiles_poly import atmload
from RadSGPoly.cooling_poly import cooling_global as cg, cooling_local as cl,\
     critical
from RadSGPoly.atmseries_poly import mcrit_disktime as mcd

def t_vs_Mc_fixed_a(delad, Y, a):

    """
    Uses already generated atmosphere profiles to calculate the atmosphere
    evolution time as a function of core mass at fixed distance in the disk.

    Input
    -----
    delad:
        adiabtic index
    Y:
        helium mass fraction
    a:
        distance in the disk

    Output
    ------
    Mc:
        core masses for the (already generated) atmospheres, in Earth masses
    t:
        cumulative cooling time ( = time for the atmosphere to become critical),
        in years
    
    """

    prms = params(Me, Re, a, delad, Y, gamma = gammafn(delad), R = Rfn(Y), \
                  Cv = Cvfn(Y, delad), Pd = Pdisk(a, mstar, FSigma, FT), \
                  Td = Tdisk(a, FT), kappa = kdust)

    t = 0 * numpy.ndarray(shape = (10), dtype = float)

    for i in range(len(t)):

        if delad == 2./7 and Y == 0.0:
            Mc = numpy.linspace(5, 14, 10)

        elif delad == 2./7 and Y == 0.3:
            Mc = numpy.linspace(1, 10, 10)
        
        atm = atmload('Mc' + str(Mc[i]) + '.npz', prms = prms)
        model = atm[0]
        param = atm[1]
        prof = atm[2]

        temp = critical(param, prof, model)
        dt = temp[-1]
        t[i] = sum(dt)

    return Mc, t / (365 * 24 * 3600)


def t_vs_a_fixed_Mc(delad, Y, Mc):


    """
    Uses already generated atmosphere profiles to calculate the atmosphere
    evolution time for a fixed core mass as a function of distance in the disk.

    Input
    -----
    delad:
        adiabtic index
    Y:
        helium mass fraction
    Mc:
        core mass in Earth masses

    Output
    ------
    a:
        array of distances in the disk; default numpy.linspace(10, 100, 10)
    t:
        cumulative cooling time ( = time for the atmosphere to become critical),
        in years for fixed core mass Mc
    
    """

    a = numpy.linspace(10, 100, 10)
    t = 0 * numpy.ndarray(shape = (len(a)), dtype = float)

    for i in range(len(a)):

        prms = params(Me, Re, a[i], delad, Y, gamma = gammafn(delad), R = Rfn(Y), \
                      Cv = Cvfn(Y, delad), Pd = Pdisk(a, mstar, FSigma, FT), \
                      Td = Tdisk(a, FT), kappa = kdust)

        atm = atmload('Mc' + str(Mc) + '.npz', prms = prms)
        model = atm[0]
        param = atm[1]
        prof = atm[2]

        temp = critical(param, prof, model)
        dt = temp[-1]
        t[i] = sum(dt)

    return a, t / (365 * 24 * 3600)


def Mc_vs_a_fixed_t(delad, Y, disklife):
    
    """

    Uses already generated atmosphere profiles to calculate the atmosphere
    evolution time for a fixed core mass as a function of distance in the disk,
    then interpolate to find the critical core mass as function of the
    distance in the disk for a fixed life time of the disk.

    Input
    -----
    delad:
        adiabatic index
    Y:
        helium mass fraction
    disklife:
        disk lifetime in years


    Output
    ------
    a:
        array of distances in the disk; default numpy.linspace(10, 100, 10)
    Mc:
        array or critical core masses for fixed disk life
        
    """

    a = numpy.linspace(10, 100, 10)
    Mcrit = 0 * numpy.ndarray(shape = (len(a)), dtype = float)
    
    for i in range(len(a)):

        Mct = t_vs_Mc_fixed_a(delad, Y, a[i])
        Mcrit[i] = mcd(Mct[0], Mct[1], disklife = disklife)

    return a, Mcrit
        
        
























