"""
This module contains functions that make generating plots easier. Used by
scripts & ipython notebooks in python/Atmscrips.
"""


import numpy
from utils.constants import Me, Re, gammafn, Rfn, Cvfn, Pdisk, Tdisk, kdust, \
     paramsEOS
from utils.parameters import FSigma, FT, mstar
import matplotlib.pyplot as plt

from RadSGRealGas.profiles import atmload
from RadSGRealGas.cooling import cooling_global as cg, cooling_local as cl,\
     critical
from RadSGRealGas.atmseries import mcrit_disktime as mcd

def t_vs_Mc_fixed_a(Y, a, returnt = 0):

    """
    Uses already generated atmosphere profiles to calculate the atmosphere
    evolution time as a function of core mass at fixed distance in the disk.

    Input
    -----
    Y:
        helium mass fraction
    a:
        distance in the disk
    returnt:
        option to also return dt  and cumuluative cooling time at each mass
        increase for every atmosphere series; set to 0 by default

    Output
    ------
    Mc:
        core masses for the (already generated) atmospheres, in Earth masses
    t:
        cumulative cooling time ( = time for the atmosphere to become critical),
        in years
    Mtotcrit:
        optional; total mass at the critical point
    dtarr:
        optional; array of instantaneous dt for each atmosphere series, in
        years
    tarr:
        optional; array of cumulative cooling times for each atmosphere series
        in years
    
    """

    prms = paramsEOS(Me, Re, Y, a, Pd = Pdisk(a, mstar, FSigma, FT), \
                  Td = Tdisk(a, FT), kappa = kdust)

    #if a == 1.0 or a == 5.0 or a  == 3.0 or a == 2.0 or a == 1.5 or a == 2.5 or a == 1.75:

    if a == 10.0:
        
        t = 0 * numpy.ndarray(shape = (11), dtype = float)
        dtarr = 0 * numpy.ndarray(shape = (11, 199), dtype = float)
        tarr = 0 * numpy.ndarray(shape = (11, 199), dtype = float)
        Mtotcrit = 0 * numpy.ndarray(shape = (199), dtype = float)

    else:

        t = 0 * numpy.ndarray(shape = (5), dtype = float)
        dtarr = 0 * numpy.ndarray(shape = (5, 199), dtype = float)
        tarr = 0 * numpy.ndarray(shape = (5, 199), dtype = float)
        Mtotcrit = 0 * numpy.ndarray(shape = (199), dtype = float)
    
    for i in range(len(t)):

        if a == 5.0:
            Mc = numpy.linspace(18, 22, 5)

        if a == 10.0:
            Mc = numpy.linspace(6, 16, 11)

	elif a == 20.0:
	    Mc = numpy.linspace(11, 15, 5)
	    
	elif a == 30.0:
	    Mc = numpy.linspace(9, 13, 5)

	elif a == 40.0:
            Mc = numpy.linspace(7, 11, 5)

	elif a == 50.0:
            Mc = numpy.linspace(6, 10, 5)
            
	elif a == 60.0 or a == 70.0 or a == 80.0:
            Mc = numpy.linspace(5, 9, 5)

        elif a == 90.0 or a == 100.0:
            Mc = numpy.linspace(4, 8, 5)

        
        atm = atmload('Mc' + str(Mc[i]) + '.npz', prms = prms)
        model = atm[0]
        param = atm[1]
        prof = atm[2]

        temp = critical(param, prof, model)
        Mtotcrit[i] = temp[0].MB[-1]
        dt = temp[-1]
        t[i] = sum(dt)

        for j in range(len(dt)):
            dtarr[i, j] = dt[j]
            tarr[i, j] = sum(dt[:j])

    if returnt == 0:
        return Mc, t / (365 * 24 * 3600)
    else:
        return Mc, t / (365 * 24 * 3600), Mtotcrit, dtarr / (365 * 24 * 3600), \
               tarr / (365 * 24 * 3600)


def t_vs_a_fixed_Mc(Y, Mc):


    """
    Uses already generated atmosphere profiles to calculate the atmosphere
    evolution time for a fixed core mass as a function of distance in the disk.

    Input
    -----
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

    a = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
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


def Mc_vs_a_fixed_t(Y, disklife):
    
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
    
    a = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    #a = [1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    Mcrit = 0 * numpy.ndarray(shape = (len(a)), dtype = float)
    
    for i in range(len(a)):

        Mct = t_vs_Mc_fixed_a(Y, a[i])
        Mcrit[i] = mcd(Mct[0], Mct[1], disklife = disklife)

    return a, Mcrit
        
        
























