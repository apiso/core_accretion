"""
This module contains functions that make generating plots easier. Used by
scripts & ipython notebooks in python/Atmscrips.
"""


import numpy
from utils.constants import Me, Re, gammafn, Rfn, Cvfn, Pdisk, Tdisk, kdust, \
     params, kdust10, kdust100
from utils.parameters import FSigma, FT, mstar
import matplotlib.pyplot as plt

from RadSGPoly.profiles_poly import atmload
from RadSGPoly.cooling_poly import cooling_global as cg, cooling_local as cl,\
     critical
from RadSGPoly.atmseries_poly import mcrit_disktime as mcd

def t_vs_Mc_fixed_a(delad, Y, a, returnt = 0, opacity = kdust):

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

    prms = params(Me, Re, a, delad, Y, gamma = gammafn(delad), R = Rfn(Y), \
                  Cv = Cvfn(Y, delad), Pd = Pdisk(a, mstar, FSigma, FT), \
                  Td = Tdisk(a, FT), kappa = kdust)

    #if prms.kappa == kdust:

    if a == 1.0 or a == 5.0 or a  == 3.0 or a == 2.0 or a == 1.5 or a == 2.5 or a == 1.75:
        t = 0 * numpy.ndarray(shape = (5), dtype = float)
        dtarr = 0 * numpy.ndarray(shape = (5, 299), dtype = float)
        tarr = 0 * numpy.ndarray(shape = (5, 299), dtype = float)
        Mtotcrit = 0 * numpy.ndarray(shape = (299), dtype = float)
    
    else:
        t = 0 * numpy.ndarray(shape = (10), dtype = float)
        dtarr = 0 * numpy.ndarray(shape = (10, 299), dtype = float)
        tarr = 0 * numpy.ndarray(shape = (10, 299), dtype = float)
        Mtotcrit = 0 * numpy.ndarray(shape = (299), dtype = float)
    
    for i in range(len(t)):

        if a == 1.0:
            if delad == 2./7 and Y == 0.0:
                Mc = numpy.linspace(16, 20, 5)
            elif delad == 2./7 and Y == 0.3 and opacity == kdust:
                Mc = numpy.linspace(10, 14, 5)
            elif delad == 2./7 and Y == 0.3 and opacity == kdust10:
                Mc = numpy.linspace(4, 8, 5)
            elif delad == 2./5 and Y == 0.3:
                Mc = numpy.linspace(40, 44, 5)

        elif a == 1.5:
            if delad == 2./5 and Y == 0.3:
                Mc = [20.0, 21.0, 22.0, 23.0, 24.0]

        elif a == 1.75:
            if delad == 2./5 and Y == 0.3:
                Mc = [15.0, 16.0, 17.0, 18.0, 19.0]
        
        elif a == 2.0 or a == 2.5:
            if delad == 2./5 and Y == 0.3:
                Mc = [8.0, 9.0, 10.0, 11.0, 12.0]

        elif a == 3.0:
            if delad == 2./5 and Y == 0.3:
                Mc = [8.0, 9.0, 10.0, 11.0, 12.0] #numpy.linspace(11, 15, 5)

        elif a == 5.0:
            if delad == 2./7 and Y == 0.0:
                Mc = numpy.linspace(11, 15, 5)
            elif delad == 2./7 and Y == 0.3 and opacity == kdust:
                Mc = numpy.linspace(7, 11, 5)
            elif delad == 2./7 and Y == 0.3 and opacity == kdust10:
                Mc = numpy.linspace(2, 6, 5)                
            elif delad == 2./5 and Y == 0.3:
                Mc = numpy.linspace(6, 10, 5)

        elif a == 10.0 and delad == 2./7 and Y == 0.3 and opacity == kdust:
                Mc = numpy.linspace(5, 14, 10)

        else:
            if delad == 2./7 and Y == 0.0:
                Mc = numpy.linspace(5, 14, 10)

            elif delad == 2./7 and Y == 0.3:
                Mc = numpy.linspace(1, 10, 10)
            
            elif delad == 2./5 and Y == 0.3:
                Mc = numpy.linspace(1, 10, 10)
    
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


##    elif prms.kappa == kdust10:
##        
##        if a == 10.0:
##            t = 0 * numpy.ndarray(shape = (5), dtype = float)
##            dtarr = 0 * numpy.ndarray(shape = (5, 299), dtype = float)
##            tarr = 0 * numpy.ndarray(shape = (5, 299), dtype = float)
##            Mtotcrit = 0 * numpy.ndarray(shape = (299), dtype = float)
##            
##        for i in range(len(t)):
##            
##            if delad == 2./7 and Y == 0.3:
##                if a == 10.0:
##                    Mc = numpy.linspace(3, 7, 5)
##
##            atm = atmload('Mc' + str(Mc[i]) + '.npz', prms = prms)
##            model = atm[0]
##            param = atm[1]
##            prof = atm[2]
##
##            temp = critical(param, prof, model)
##            Mtotcrit[i] = temp[0].MB[-1]
##            dt = temp[-1]
##            t[i] = sum(dt)
##
##            for j in range(len(dt)):
##                dtarr[i, j] = dt[j]
##                tarr[i, j] = sum(dt[:j])


    if opacity == kdust10:
        t = t / 10
    elif opacity == kdust100:
        t = t / 100
        
    if returnt == 0:
        return Mc, t / (365 * 24 * 3600)
    else:
        return Mc, t / (365 * 24 * 3600), Mtotcrit, dtarr / (365 * 24 * 3600), \
               tarr / (365 * 24 * 3600)


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


def Mc_vs_a_fixed_t(delad, Y, disklife, opacity = kdust):
    
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
    if delad == 2./5:
    	a = [1.0, 1.5, 1.75, 2.0, 2.5, 3.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0] # numpy.linspace(10, 100, 10)
    else:
	a = [1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    Mcrit = 0 * numpy.ndarray(shape = (len(a)), dtype = float)
    
    for i in range(len(a)):

        Mct = t_vs_Mc_fixed_a(delad, Y, a[i], opacity = opacity)
        Mcrit[i] = mcd(Mct[0], Mct[1], disklife = disklife)

    return a, Mcrit
        
        
























