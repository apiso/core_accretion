"""

Generates series of atmosphere profiles for a given core mass and semi-major axis, and writes them in .npz files.
Also contains functions to load profiles.

"""



from utils.constants import G, kb, mp, Rb, Me, Re, RH, RHe, RHill, sigma
from utils.parameters import rc, Mc, Y, beta, R, gamma, delad, Cv, a
from utils.userpath import userpath
import numpy
import scipy
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate
from utils.table_size import Nt, Np, Ncol
from utils.disk_model import Td, Pd, Sd, rhod, kappa
from utils import cb_corrections as cb
from utils.interpolation_functions_logp import interpolations
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.optimize import fminbound, brentq
#from shooting_poly import Tcore, shoot, minent
from shooting_poly_Pc import Pcore, shoot, minent
from utils.myspace import loglinspace


def profiles_write(n, nent, T1, T2, log10S1, log10S2, paramfilename, proffilename, tol = 10**(-25), savefile = 1, partial = 0):

    if partial == 0:
        log10Smin = minent(log10S1, log10S2, T1, T2, n, tol = 10**(-25), tols = 10**(-5))
        log10Smax = Sd - 10**(-5) * Sd
    else:
        log10Smin = log10S1
        log10Smax = log10S2
    
    prof = numpy.recarray(shape = (2 * nent, n), dtype = [('t', float), ('r', float), ('m', float), ('P', float), ('rho', float)])
    param = numpy.recarray(shape = (2 * nent), dtype = [('log10S', float), ('Mcb', float), ('Mrad', float), ('rcb', float), ('RB', float), \
                                                     ('RHill', float), ('Pc', float), ('Pcb', float), ('Tc', float), \
                                                     ('Tcb', float), ('Egcb', float), ('Ucb', float), ('Etotcb', float), \
                                                     ('EgHill', float), ('UHill', float), ('EtotHill', float),
                                                     ('L', float), ('vircheck', float), ('err', float)])
    
    log10Shigh = numpy.linspace(log10Smin, log10Smax, nent)
    log10Slow = log10Shigh[::-1]

    for i in range(nent):
        sol = shoot(log10Shigh[i], T1, T2, n, tol)
        sollow = sol[0][0]
        solhigh = sol[1][0]
        errlow = sol[0][1]
        errhigh = sol[1][1]
        jlow = nent - i -1
        jhigh = i + nent
        
        param.Mcb[jlow], param.Mrad[jlow], param.rcb[jlow], param.RB[jlow], param.RHill[jlow], param.Pc[jlow], \
                      param.Pcb[jlow], param.Tc[jlow], param.Tcb[jlow], param.Egcb[jlow], param.Ucb[jlow], param.Etotcb[jlow], \
                      param.EgHill[jlow], param.UHill[jlow], param.EtotHill[jlow], \
                      param.L[jlow], param.vircheck[jlow] = sollow[5:]
        param.log10S[jlow] = log10Slow[jlow]
        param.err[jlow] = errlow
        
        param.Mcb[jhigh], param.Mrad[jhigh], param.rcb[jhigh], param.RB[jhigh],\
                    param.RHill[jhigh], param.Pc[jhigh], param.Pcb[jhigh], param.Tc[jhigh], \
                    param.Tcb[jhigh], param.Egcb[jhigh], param.Ucb[jhigh], param.Etotcb[jhigh], \
                    param.EgHill[jhigh], param.UHill[jhigh], param.EtotHill[jhigh], \
                      param.L[jhigh], param.vircheck[jhigh] = solhigh[5:]
        param.log10S[jhigh] = log10Shigh[i]
        param.err[jhigh] = errhigh

        for k in range(n):
            prof.t[jlow, k], prof.r[jlow, k], prof.m[jlow, k], prof.P[jlow, k], prof.rho[jlow, k] = \
                         sollow[0][k], sollow[1][k], sollow[2][k], sollow[3][k], sollow[4][k]
            prof.t[jhigh, k], prof.r[jhigh, k], prof.m[jhigh, k], prof.P[jhigh, k], prof.rho[jhigh, k] = \
                         solhigh[0][k], solhigh[1][k], solhigh[2][k], solhigh[3][k], solhigh[4][k]

    if savefile == 1:
        paramfilename = userpath + '/dat/MODELS/RadNoSGPoly/' + paramfilename
        numpy.savez_compressed(paramfilename, param = param)

        proffilename = userpath + '/dat/MODELS/RadNoSGPoly/' + proffilename
        numpy.savez_compressed(proffilename, prof = prof)
    
    return param, prof


def paramload(filename):
    npzdat = numpy.load(userpath + '/dat/MODELS/RadNoSGPoly/' + filename)
    param = npzdat['param'].view(numpy.recarray)

    return param

def profload(filename):
    npzdat = numpy.load(userpath + '/dat/MODELS/RadNoSGPoly/' + filename)
    prof = npzdat['prof'].view(numpy.recarray)

    return prof
