"""

Computes global and local energy balance for a given series of atmospheres.

"""


from utils.constants import G, kb, mp, Rb, Me, Re, RH, RHe, sigma
from utils.parameters import rc, Mc, Y, beta, R, gamma, delad, Cv, a
import numpy
import scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate
from utils.table_size import Nt, Np, Ncol
from utils.disk_model import Td, Pd, rhod, kappa
from utils import cb_corrections as cb
from scipy.integrate import odeint
from scipy.interpolate import interp1d

#---------------------------------------------------------------------------


def cooling_global(atmset, atmprofile):
    
    """
    Determines the time evolution of a given atmosphere. Cooling time is calculated
    based on the following formula:
    dt = - (dE - <ecb> dm + <Pcb> dV + <G*M/(3*R)> dm), where dV is the volume change at constant mass

    Input:
        -atmset: set of atmosphere parameters; read from dat file of atmosphere parameters with
        atmparamread function from wr_selfgrav_poly module
        -atmprofile: set of atmosphere profiles; read from dat file of atmosphere profiles with
        atmprofileread function from wr_selfgrav_poly module
        
    Returns:
        -dt for each consecutive sets of atmospheres
        -total time (sum(dt)) in units of 3 * 10^6 years (typical disk lifetime)
    """
    
    L = atmset.L
    deltav = 0 * numpy.ndarray(shape = (len(L) - 1), dtype = float)

    Tcb = atmset.Tcb
    Pcb = atmset.Pcb #array of pressures at the conv. boundary for each atm.
    rcb = atmset.rcb * Re #array of radii of the conv. boundary
    Mcb = atmset.Mcb * Me #array of masses of the conv. region 
    Tc = atmset.Tc #array of core temperatures
    Pc = atmset.Pc #array of core pressures
    Etot = atmset.Etotcb #array of total energies in the conv. region
    ecb = Cv * Tcb - G * Mcb / rcb
            #ecb = specific energy of accreted mass; ecb = ucb - G * M / rb,
            #with ucb the internal energy per unit mass at the convective boundary
    deltamcb = Mcb[1:] - Mcb[:-1]
    deltae = Etot[1:] - Etot[:-1]
    ecbav = (ecb[1:] + ecb[:-1]) / 2
    Lav = (L[1:] + L[:-1]) / 2
    Pbav = (Pcb[1:] + Pcb[:-1]) / 2
    
    for i in range(len(L) - 1):
        
        if (deltamcb[i] > 0):
            
            mnewprof = atmprofile.m[i+1] #mass profile of 'new' atmosphere (i+1)
            rnewprof = atmprofile.r[i+1]#radius profile of 'new' atmosphere (i+11)
            f = interp1d(mnewprof, rnewprof) #interpolation function
            rcbold = f(Mcb[i]) #we find the radius in the new atmosphere (i+1) at which
                              #the mass is equal to the mass of the old atmosphere (i)
                              #in order to be able to calculate dV at constant m
            deltav[i] = (4 * numpy.pi / 3) * (rcbold**3 - rcb[i]**3) #dV
            
        else:
            i = i + 1

    t = ( - deltae + ecbav * deltamcb - Pbav * deltav) \
        / Lav
    #tid = - deltae / Lav #time evolution if we only take into account dE
    return t, sum(t / (365 * 24 * 3600)) / (3 * 10**6), \
           deltae, ecbav * deltamcb, Pbav * deltav


def cooling_local(param, prof):

    t = cooling_global(param, prof)[0]
    n = numpy.shape(param)[0]
    npoints = numpy.shape(prof)[1]
    M = param.Mcb * Me
    Mav = (M[1:] + M[:-1]) / 2
    Ldt = 0 * numpy.ndarray(shape = (n - 1), dtype = float)
    
    for i in range(n - 1):

        if M[i + 1] > M[i]:
        
            m = numpy.linspace(Mc, M[i], npoints)
            mmid = (m[:-1] + m[1:]) / 2
            dm = m[1] - m[0]
            
            m1 = prof.m[i]
            T1 = prof.t[i]
            P1 = prof.P[i]

            m2 = prof.m[i + 1]
            T2 = prof.t[i + 1]
            P2 = prof.P[i + 1]

            fP1 = interp1d(m1, P1)
            fT1 = interp1d(m1, T1)
            
            fP2 = interp1d(m2, P2)
            fT2 = interp1d(m2, T2)

            P1int = fP1(mmid)
            P2int = fP2(mmid)
            T1int = fT1(mmid)
            T2int = fT2(mmid)

            Tav = (T1int + T2int) / 2

            dS = (R / delad) * numpy.log((T2int / T1int) * (P1int / P2int)**delad)
            Ldt[i] = - sum(Tav * dS * dm)
            
        else:
            i = i + 1

    return Ldt









