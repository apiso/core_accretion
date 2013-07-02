"""

Computes global and local energy balance for a given series of atmospheres.
Also used to find the critical core mass.

"""


from utils.constants import G, kb, mp, Rb, Me, Re, Msun, RH, RHe, sigma, \
     cmperau, RHill, gammafn, mufn, Rfn, Cvfn, kdust, Tdisk, Pdisk, params
import numpy
import scipy
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import namedtuple
from scipy import integrate
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from shooting_poly import Ltop, shoot, prms

#---------------------------------------------------------------------------

def cooling_global(atmset, atmprofile, prms = prms, out = 'rcb'):
    
    """
    Determines the time evolution of a given atmosphere. Cooling time is
    calculated based on the following formula:
    dt = - (dE - <ecb> dm + <Pcb> dV  dm) / L,
    where dV is the volume change at constant mass

    Input
    -----
    atmset:
        the 'param' recarray (see profiles_poly.py)
    atmprofile:
        the 'prof' recarray (see profiles_poly.py)
    prms:
        gas, disk and core parameters; default prms imported from shoot
    out:
        specifies the surfaces where cooling is computed; default RCB
        

    Output
    ------
    t:
        array of delta t's between subsequent atmospheres in seconds
    t_cumulative:
        cumulative cooling time in units of 3 Myrs
    deltae, eaccav * deltamout, Pout * deltav:
        the cooling terms in the energy equation
        
    """
    
    n = len(atmset)
    L = atmset.L
    deltav = 0 * numpy.ndarray(shape = (len(L) - 1), dtype = float)

    if out == 'rcb':
        Tout = atmset.Tcb
        Pout = atmset.Pcb 
        rout = atmset.rcb * Re
        Mout = atmset.Mcb * Me 
        Eout = atmset.Etotcb 

    elif out == 'RHill':
        Tout = numpy.array([prms.Td] * n)
        Pout = numpy.array([prms.Pd] * n) 
        rout = atmset.RHill * Re 
        Mout = atmset.Mtot * Me 
        Eout = atmset.EtotHill 

    elif out == 'RB':
        Tout = atmset.TB
        Pout = atmset.PB 
        rout = atmset.RB * Re 
        Mout = atmset.MB * Me 
        Eout = atmset.EtotB 
        
    eaccout = prms.Cv * Tout - G * Mout / rout
    deltamout = Mout[1:] - Mout[:-1]
    deltae = Eout[1:] - Eout[:-1]
    eaccav = (eaccout[1:] + eaccout[:-1]) / 2
    Lav = (L[1:] + L[:-1]) / 2
    Poutav = (Pout[1:] + Pout[:-1]) / 2
    
    
    for i in range(len(L) - 1):
        
        if (deltamout[i] > 0):
            
            mnewprof = atmprofile.m[i+1] #mass profile of 'new' atmosphere (i+1)
            rnewprof = atmprofile.r[i+1]#radius profile of 'new' atmosphere(i+1)
            f = interp1d(mnewprof, rnewprof) #interpolation function
            routold = f(Mout[i]) #we find the radius in the new atmosphere (i+1)
                                 #at which the mass is equal to the mass of the
                                 #old atmosphere (i) in order to be able to
                                 #calculate dV at constant m
            
            deltav[i] = (4 * numpy.pi / 3) * (routold**3 - rout[i]**3)
            
        else:
            i = i + 1

    t = ( - deltae + eaccav * deltamout - Poutav * deltav) \
            / Lav
    
    return t, sum(t / (365 * 24 * 3600)) / (3 * 10**6), \
          deltae, eaccav * deltamout, Poutav * deltav



def cooling_local(param, prof, prms = prms, out = 'rcb', onlyrad = 0):

    """
    Calculates cooling terms from dL/dm = -T dS/dt


    Input
    -----
    atmset:
        the 'param' recarray (see profiles_poly.py)
    atmprofile:
        the 'prof' recarray (see profiles_poly.py)
    prms:
        gas, disk and core parameters; default prms imported from shoot
    out:
        specifies the surfaces where cooling is computed; default RCB
    onlyrad:
        flag to only calculate the luminosity in the radiative region
        

    Output
    ------
    Ldt:
        Ldt = integral(-T dS dm)

    """
    
    n = numpy.shape(param)[0]
    npoints = numpy.shape(prof)[1]

    Mcb = param.Mcb * Me
    Mcbav = (Mcb[1:] + Mcb[:-1]) / 2
    Mtot = param.Mtot * Me
    MB = param.MB * Me
    
    if out == 'rcb':
        M = Mcb
    elif out == 'RHill':
        M = Mtot
    elif out == 'RB':
        M = MB
        
    Mav = (M[1:] + M[:-1]) / 2
    Ldt = 0 * numpy.ndarray(shape = (n - 1), dtype = float)
    
    for i in range(n - 1):

        if M[i + 1] > M[i]:

            if onlyrad == 0:
                m = numpy.linspace(prms.Mco, M[i], npoints)
            elif onlyrad != 0 and out == 'RHill':
                m = numpy.linspace(Mcbav[i], Mtot[i], npoints)
            elif onlyrad != 0 and out == 'RB':
                m = numpy.linspace(Mcbav[i], MB[i], npoints)
            else:
                print "Wrong choice of boundary."
                sys.exit()

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

            success = 0
        
            while success != 1:
                try:
                    P1int = fP1(mmid)
                    P2int = fP2(mmid)
                    T1int = fT1(mmid)
                    T2int = fT2(mmid)

                    success = 1

                except ValueError:
                    mmid = mmid[1:]

            Tav = (T1int + T2int) / 2

            dS = (prms.R / prms.delad) * numpy.log((T2int / T1int) * \
                                                   (P1int / P2int)**prms.delad)
            Ldt[i] = - sum(Tav * dS * dm)
            
        else:
            i = i + 1

    return Ldt

def critical(param, prof, prms = prms):
    
    """

    For a given atmosphere profiles, finds where the atmosphere becomes
    critical. This is defined by the minimum between mass doubling and
    entropy minimum
    
    Input
    -----
    atmset:
        the 'param' recarray (see profiles_poly.py)
    atmprofile:
        the 'prof' recarray (see profiles_poly.py)
    prms:
        gas, disk and core parameters; default prms imported from shoot


    Output
    ------
    param, prof:
        the atmosphere recarrays up to the critical point
    t:
        the delta t between two subsequent atmospheres up to the critical point
    
    """
    
    dt = cooling_global(param, prof, prms)[0]
    Lav = (param.L[:-1] + param.L[1:]) / 2
    
    for i in range(len(dt)):
        if dt[i] < 0:
            break
    for j in range(len(param.MB)):
        if param.MB[j] > 2 * prms.Mco / Me:
            break
        
    k = numpy.array([i, j]).min()
    return param[:k], prof[:k], dt[:k]




def tacc(param, prof, prms):

    x = critical(param, prof, prms)
    param, prof, dttrue = x
    L = param.L * 100

    l = list(L)
    Lacc = L[1 : l.index(L.min()) + 1]
    #Lacc = numpy.logspace(numpy.log10(L[0] / 5), \
    #                             numpy.log10(L.min()), 100)
    dt = 0 * numpy.ndarray(shape = (len(Lacc), len(L) - 1), dtype = float)
    t = 0 * numpy.ndarray(shape = (len(Lacc)), dtype = float)
    Macc = 0 * numpy.ndarray(shape = (len(Lacc)), dtype = float)

    #l = list(L)
    #Lnew = L[:l.index(L.min()) + 1]
    #f = interp1d(Lnew[::-1], Lacc[::-1])
    #Laccnew = f(Lnew)
        
    for i in range(len(Lacc)):
        for j in range(len(L) - 1):
        
            dt[i, j] = L[j]/ (L[j] - Lacc[i]) * dttrue[j] / 100
            for k in range(len(L) - 1):
                if math.isinf(dt[i, k]):
                    break
            t[i] = sum(dt[i, :k]) / (365 * 24 * 3600 * 10**6)
            #for k in range(k, len(L) - 1):
        dt[i, k:] = 0
        Macc[i] = param.MB[k]
    return dt, param.MB, t, Macc, Lacc

def tacc_plot(param, prof, prms, n):
    x = tacc(param, prof, prms)
    dt = x[0] / (365 * 24 * 3600)
    MB = x[1]

    crit = critical(param, prof, prms)
    paramcrit = crit[0]
    dttrue = crit[2]
    ttrue = 0 * numpy.ndarray(shape = len(dttrue), dtype = float)
    for i in range(len(ttrue)):
        ttrue[i] = sum(dttrue[:i]) / (365 * 24 * 3600) / 100

    tcum = 0 * numpy.ndarray(shape = (numpy.shape(dt)[0], numpy.shape(dt)[1]), \
                       dtype = float)
    for i in range(numpy.shape(tcum)[0]):
        for j in range(numpy.shape(tcum)[1]):
            tcum[i, j] = sum(dt[i, :j])

    for i in range(numpy.shape(tcum)[0]):
        for j in range(numpy.shape(tcum)[1] - 1):
            if tcum[i, j - 1] == tcum[i, j]:
                break
        tcum[i, j:] = 0
            
    
    rnb = numpy.linspace(0, 256, numpy.shape(dt)[0] / n)
    
    for i in range(numpy.shape(dt)[0] / n):
        plt.semilogy(MB[:-1], tcum[n * i], c = cm.rainbow(numpy.int(rnb[i])))
    plt.semilogy(paramcrit.MB, ttrue, '--', color = 'black')
    plt.xlabel(r'$M_{\rm{atm}}\,[M_{\oplus}]$')
    plt.ylabel(r'$t_{\rm{elapsed}}$ [yrs]')
    plt.xlim(xmax = 3.8)
    plt.show()





