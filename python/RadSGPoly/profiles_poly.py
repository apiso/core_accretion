"""

Generates series of atmosphere profiles for a set of gas, disk and core
conditions and writes them in .npz files. Also contains functions to load
profiles.

"""



from utils.constants import G, kb, mp, Rb, Me, Re, Msun, RH, RHe, sigma, \
     cmperau, RHill, gammafn, mufn, Rfn, Cvfn, kdust, Tdisk, Pdisk, params, \
     kdust, kdustbeta1, kdust10, kdust100
from utils.userpath import userpath
import numpy
import scipy
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from types import FunctionType as function
from collections import namedtuple
from scipy import integrate
from utils import cb_corrections as cb
from utils.interpolation_functions_logp import interpolations
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.optimize import fminbound, brentq
from shooting_poly import Ltop, shoot, prms


def profiles_write(n, nMpoints, L1, L2, Mmin, Mmax, filename, prms = prms, \
                   tol = 10**(-25), log = 1, savefile = 1, disk = 1):

    """
    Uses the output of shoot (described in shooting_poly) to generate recarrays
    and write .npz atmosphere files.

    Input
    -----
    n:
        number of grid points for each atmosphere of fixed mass
    nMpoints:
        number of atmosphere profiles of fixed mass for the given core mass
    L1, L2:
        trial luminosity values. L1 and L2 must bracket the solution.
    Mmin, Mmax:
        minimum and maximum atmosphere + core masses in g
    filename:
        name of the file in which the data is saved; typically 'Mcx',
        where x is the core mass in Earth masses
    prms:
        gas, disk and core parameters; default prms imported from shoot
    tol:
        relative luminosity matching error; default 10**(-25)
    log:
        flag used to specify whether the mass spacing should be logarithmic
        or linear; default 1 (spacing is logarithmic)
    savefile:
        flag used to specify whether the recarray should be saved to a file
        or not; default 1, file is saved


    Output
    ------
    model:
        array of gas, disk and core parameters
    param:
        recarray of
            Mtot:
                total atmosphere mass in Earth masses
            Mcb:
                total mass at the radiative convective boundary in Earth masses
            MB:
                total mass inside the Bondi radius in Earth masses; this is
                what we define as the mass of the atmosphere
            rcb, RB, RHill:
                radiative convective boundary, Bondi and Hill radii in Earth
                radii
            Pc, Pcb, PB:
                pressure at the core, RCB and RB
            Tc, Tcb, TB:
                temperature at the core, RCB and RB
            Egcb, Ucb, Etotcb:
                gravitational, internal and total energy at the RCB in ergs
            EgB, UB, EtotB:
                gravitational, internal and total energy at the Bondi radius
                in ergs
            EgHill, UHill, EtotHill:
                gravitational, internal and total energy at the Hill radius
                in ergs
            L:
                matched luminosity in ergs/s
            vircb, virHill:
                checks of virial equilibrium at RCB and RHill
            err:
                luminosity matching error
    prof:
        recarray of r, P(r), T(r), m(r), rho(r), delrad(r)
        
    
    """

    model = numpy.array([(prms.Mco, prms.rco, prms.a, prms.delad, prms.Y, \
                          prms.gamma, prms.R, prms.Cv, prms.Pd, prms.Td, \
                          prms.kappa)], \
                        dtype = [('Mco', float), ('rco', float), ('a', float), \
                                 ('delad', float), ('Y', float), \
                                 ('gamma', float), \
                                 ('R', float), ('Cv', float), ('Pd', float), \
                                 ('Td', float), ('kappa', function)])
    
    prof = numpy.recarray(shape = (nMpoints, n), \
                          dtype = [('r', float), ('P', float), ('t', float), \
                                   ('m', float), ('rho', float), ('delrad', float)])
    param = numpy.recarray(\
        shape = (nMpoints), \
        dtype = [('Mtot', float), ('Mcb', float), ('MB', float), ('rcb', float),\
                 ('RB', float), ('RHill', float), ('Pc', float), ('Pcb', float),\
                 ('PB', float), ('Tc', float), ('Tcb', float), ('TB', float), \
                 ('Egcb', float), ('Ucb', float), ('Etotcb', float), \
                 ('EgB', float), ('UB', float), ('EtotB', float), \
                 ('EgHill', float), ('UHill', float), ('EtotHill', float), \
                 ('L', float), ('vircb', float), ('virHill', float), \
                 ('err', float)])
        
    if log == 1:
        mass = numpy.logspace(numpy.log10(Mmin), numpy.log10(Mmax), nMpoints)
    else:
        mass = numpy.linspace(Mmin, Mmax, nMpoints)

    for i in range(nMpoints):
        sol = shoot(mass[i], L1, L2, n, tol, prms)
        
        param.Mtot[i], param.Mcb[i], param.MB[i], param.rcb[i], param.RB[i], \
                       param.RHill[i], param.Pc[i], param.Pcb[i], param.PB[i],\
                       param.Tc[i], param.Tcb[i], param.TB[i], param.Egcb[i], \
                       param.Ucb[i], param.Etotcb[i], param.EgB[i], \
                       param.UB[i], param.EtotB[i], param.EgHill[i], \
                       param.UHill[i], param.EtotHill[i], param.L[i], \
                       param.vircb[i], param.virHill[i], param.err[i] = sol[6:]

        for k in range(n):
            prof.r[i, k], prof.P[i, k], prof.t[i, k], prof.m[i, k], \
                      prof.rho[i, k], prof.delrad[i, k] = \
                      sol[0][k], sol[1][k], sol[2][k], sol[3][k], sol[4][k], \
                      sol[5][k]
        print i


    if savefile == 1:
        if disk == 1:
            if prms.kappa == kdust or prms.kappa == kdust10:
                paramfilename = userpath + '/dat_ana/MODELS/RadSGPoly/' + 'delad' + \
                                str(prms.delad)[:4] + '/Y' + str(prms.Y) + '/' + \
                                str(prms.a) + 'AU/' + filename
                numpy.savez_compressed(paramfilename, model = model, param = param, \
                                   prof = prof)
##            elif prms.kappa == kdust10:
##                paramfilename = userpath + '/dat_ana/MODELS/RadSGPoly/' + 'delad' + \
##                                str(prms.delad)[:4] + '/Y' + str(prms.Y) + '/kdust10/' + \
##                                str(prms.a) + 'AU/' + filename
                numpy.savez_compressed(paramfilename, model = model, param = param, \
                                   prof = prof)
            elif prms.kappa == kdust100:
                paramfilename = userpath + '/dat_ana/MODELS/RadSGPoly/' + 'delad' + \
                                str(prms.delad)[:4] + '/Y' + str(prms.Y) + '/kdust100/' + \
                                str(prms.a) + 'AU/' + filename
                numpy.savez_compressed(paramfilename, model = model, param = param, \
                                   prof = prof)
            elif prms.kappa == kdustbeta1:
                paramfilename = userpath + '/dat_ana/MODELS/RadSGPoly/' + 'delad' + \
                                str(prms.delad)[:4] + '/Y' + str(prms.Y) + '/kdustbeta1/' + \
                                str(prms.a) + 'AU/' + filename
                numpy.savez_compressed(paramfilename, model = model, param = param, \
                                   prof = prof)

        else:
            paramfilename = userpath + '/dat_ana/MODELS/RadSGPoly/' + 'Td' + \
                            str(prms.Td)[:6] + '_Pd' + str(prms.Pd)[:7] + \
                            'Mc' + str(prms.Mco/Me)
            numpy.savez_compressed(paramfilename, model = model, param = param, \
                               prof = prof)
    
    return model, param, prof


def atmload(filename, prms = prms, disk = 1):
    if disk == 1:
        if prms.kappa == kdust or prms.kappa == kdust10:
            npzdat = numpy.load(userpath + '/dat_ana/MODELS/RadSGPoly/' + 'delad' + \
                            str(prms.delad)[:4] + '/Y' + str(prms.Y) + '/' + \
                            str(prms.a) + 'AU/' + filename)
##        elif prms.kappa == kdust10:
##            npzdat = numpy.load(userpath + '/dat_ana/MODELS/RadSGPoly/' + 'delad' + \
##                                str(prms.delad)[:4] + '/Y' + str(prms.Y) + '/kdust10/' + \
##                                str(prms.a) + 'AU/' + filename)
        elif prms.kappa == kdust100:
            npzdat = numpy.load(userpath + '/dat_ana/MODELS/RadSGPoly/' + 'delad' + \
                                str(prms.delad)[:4] + '/Y' + str(prms.Y) + '/' + \
                                str(prms.a) + 'AU/' + filename)
        elif prms.kappa == kdustbeta1:
            npzdat = numpy.load(userpath + '/dat_ana/MODELS/RadSGPoly/' + 'delad' + \
                                str(prms.delad)[:4] + '/Y' + str(prms.Y) + '/kdustbeta1/' + \
                                str(prms.a) + 'AU/' + filename)  

    else:
        npzdat = numpy.load(userpath + '/dat_ana/MODELS/RadSGPoly/' + 'Td' + \
                            str(prms.Td)[:6] + '_Pd' + str(prms.Pd)[:7] + \
                            'Mc' + str(prms.Mco/Me) + '.npz')
        
    model = npzdat['model'].view(numpy.recarray)
    param = npzdat['param'].view(numpy.recarray)
    prof = npzdat['prof'].view(numpy.recarray)

    npzdat.close()

    return model, param, prof


