"""
This module contains functions that make generating plots easier. Used by
scripts & ipython notebooks in python/Atmscrips.
"""


import numpy
from utils.constants import Me, Re, gammafn, Rfn, Cvfn, Pdisk, Tdisk, kdust, \
     paramsEOS, kdustbeta1
from utils.parameters import FSigma, FT, mstar
import matplotlib.pyplot as plt
import types
from RadSGRealGas.profiles import atmload
from RadSGRealGas.cooling import cooling_global as cg, cooling_local as cl,\
     critical, cb_options
from RadSGRealGas.atmseries import mcrit_disktime as mcd
from RadSGRealGas.gg_opacity import interp_opacity

def t_vs_Mc_fixed_a(Y, a, returnt = 0, opacity = kdust):

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
                  Td = Tdisk(a, FT), kappa = opacity)

    if prms.kappa == kdust:

    #if a == 1.0 or a == 5.0 or a  == 3.0 or a == 2.0 or a == 1.5 or a == 2.5 or a == 1.75:

        if a == 10.0:
            
            t = 0 * numpy.ndarray(shape = (10), dtype = float)
            dtarr = 0 * numpy.ndarray(shape = (10, 199), dtype = float)
            tarr = 0 * numpy.ndarray(shape = (10, 199), dtype = float)
            Mtotcrit = 0 * numpy.ndarray(shape = (199), dtype = float)
            MBondi = 0 * numpy.ndarray(shape = (10, 199), dtype = float)

        else:

            t = 0 * numpy.ndarray(shape = (5), dtype = float)
            dtarr = 0 * numpy.ndarray(shape = (5, 199), dtype = float)
            tarr = 0 * numpy.ndarray(shape = (5, 199), dtype = float)
            Mtotcrit = 0 * numpy.ndarray(shape = (199), dtype = float)
            MBondi = 0 * numpy.ndarray(shape = (5, 199), dtype = float)
        
        for i in range(len(t)):

            if a == 5.0:
                Mc = numpy.linspace(18, 22, 5)

            if a == 10.0:
                Mc = numpy.linspace(7, 16, 10)

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
            MBonditemp = temp[0].MB
            dt = temp[-1]
            t[i] = sum(dt)

            for j in range(len(dt)):
                dtarr[i, j] = dt[j]
                tarr[i, j] = sum(dt[:j])
                MBondi[i, j] = MBonditemp[j]

    elif prms.kappa == kdustbeta1:

        if a <= 30.0:
            
            t = 0 * numpy.ndarray(shape = (3), dtype = float)
            dtarr = 0 * numpy.ndarray(shape = (3, 99), dtype = float)
            tarr = 0 * numpy.ndarray(shape = (3, 99), dtype = float)
            Mtotcrit = 0 * numpy.ndarray(shape = (99), dtype = float)
            MBondi = 0 * numpy.ndarray(shape = (3, 99), dtype = float)

        else:

            t = 0 * numpy.ndarray(shape = (4), dtype = float)
            dtarr = 0 * numpy.ndarray(shape = (4, 99), dtype = float)
            tarr = 0 * numpy.ndarray(shape = (4, 99), dtype = float)
            Mtotcrit = 0 * numpy.ndarray(shape = (99), dtype = float)
            MBondi = 0 * numpy.ndarray(shape = (4, 99), dtype = float)
        
        for i in range(len(t)):

            if a <= 30.0:
                Mc = numpy.linspace(3, 5, 3)

            else:
                Mc = numpy.linspace(2, 5, 4)

            

        
            atm = atmload('Mc' + str(Mc[i]) + '.npz', prms = prms)
            model = atm[0]
            param = atm[1]
            prof = atm[2]

            temp = critical(param, prof, model)
            Mtotcrit[i] = temp[0].MB[-1]
            MBonditemp = temp[0].MB
            dt = temp[-1]
            t[i] = sum(dt)

            for j in range(len(dt)):
                dtarr[i, j] = dt[j]
                tarr[i, j] = sum(dt[:j])
                MBondi[i, j] = MBonditemp[j]

    elif prms.kappa == interp_opacity:
            
        t = 0 * numpy.ndarray(shape = (3), dtype = float)
        dtarr = 0 * numpy.ndarray(shape = (3, 199), dtype = float)
        tarr = 0 * numpy.ndarray(shape = (3, 199), dtype = float)
        Mtotcrit = 0 * numpy.ndarray(shape = (199), dtype = float)
        MBondi = 0 * numpy.ndarray(shape = (3, 199), dtype = float)

        
        for i in range(len(t)):

            Mc = numpy.linspace(2, 4, 3)            
        
            atm = atmload('Mc' + str(Mc[i]) + '.npz', prms = prms)
            model = atm[0]
            param = atm[1]
            prof = atm[2]

            temp = critical(param, prof, model)
            Mtotcrit[i] = temp[0].MB[-1]
            MBonditemp = temp[0].MB
            dt = temp[-1]
            t[i] = sum(dt)

            for j in range(len(dt)):
                dtarr[i, j] = dt[j]
                tarr[i, j] = sum(dt[:j])
                MBondi[i, j] = MBonditemp[j]
        

    if returnt == 0:
        return Mc, t / (365 * 24 * 3600)
    else:
        return Mc, t / (365 * 24 * 3600), Mtotcrit, dtarr / (365 * 24 * 3600), \
               tarr / (365 * 24 * 3600), MBondi


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


def Mc_vs_a_fixed_t(Y, disklife, opacity = kdust):
    
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
    
    a = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    #a = [1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    Mcrit = 0 * numpy.ndarray(shape = (len(a)), dtype = float)
    
    for i in range(len(a)):

        Mct = t_vs_Mc_fixed_a(Y, a[i], opacity = opacity)
        Mcrit[i] = mcd(Mct[0], Mct[1], disklife = disklife)

    return a, Mcrit
        
        
def time_evol(Y, a, Mc, p = 3.5):

    rhoc = 3.2
    rc = (3 * Mc / (4 * numpy.pi * rhoc))**(1./3)
    
    prms = paramsEOS(Mc, rc, Y, a, Pd = Pdisk(a, mstar, FSigma, FT), \
                Td = Tdisk(a, FT), kappa = interp_opacity)

    atm = atmload('Mc' + str(Mc/Me) + '.npz', prms = prms, p = p)
    model, param, prof = atm

    opt = cb_options(param) 

    if type(opt) == types.NoneType:

        t = 0 * numpy.ndarray(shape = (1), dtype = float)
        tcum = 0 * numpy.ndarray(shape = (200), dtype = float)
        MB = 0 * numpy.ndarray(shape = (200), dtype = float)

        MBf = []
        
        temp = critical(param[1:], prof[1:], model, out = 'rcb', outrad = opt)
        paramcrit, profcrit, dt = temp   
        MBf = numpy.append(MBf, paramcrit.MB[-1])
        t = sum(dt)

        for j in range(len(dt)):
            tcum[j] = sum(dt[:j])
            MB[j] = paramcrit.MB[j]

        
    else:
        
        t = 0 * numpy.ndarray(shape = (len(opt)), dtype = float)
        tcum = 0 * numpy.ndarray(shape = (len(opt), 200), dtype = float)
        MB = 0 * numpy.ndarray(shape = (len(opt), 200), dtype = float)

        MBf = []
        
        for i in range(len(opt)):

            if a != 10.0 and Mc != 4*Me:
                temp = critical(param, prof, model, out = 'rcb', outrad = opt[i])
            else:
                temp = critical(param[1:], prof[1:], model, out = 'rcb', outrad = opt[i])
            paramcrit, profcrit, dt = temp   
            MBf = numpy.append(MBf, paramcrit.MB[-1])
            t[i] = sum(dt)

            for j in range(len(dt)):
                tcum[i, j] = sum(dt[:j])
                MB[i, j] = paramcrit.MB[j]

    tRB = []
    
    dtRB = cg(paramcrit, profcrit, model, out = 'RB')[0]
    for k in range(len(dtRB)):
        tRB = numpy.append(tRB, sum(dtRB[:k]))
    
    #dtRH = cg(paramcrit, profcrit, model, out = 'RHill')[0]
    #for k in range(len(dtRH)):
    #    tRH = numpy.append(tRH, sum(dtRH[:k]))                                


    return paramcrit.MB, tRB, MB, tcum, MBf, t, param, opt





















