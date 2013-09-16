"""

Computes global and local energy balance for a given series of atmospheres.
Also used to find the critical core mass.

"""


from utils.constants import G, kb, mp, Rb, Me, Re, Msun, RH, RHe, sigma, \
     cmperau, RHill, gammafn, mufn, Rfn, Cvfn, kdust, Tdisk, Pdisk, params
from utils.interpolation_functions import interplog10u, interplog10S
import numpy
import scipy
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import namedtuple
from scipy import integrate
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from shooting import prms
from RadSGRealGas.gg_opacity import interp_opacity

#---------------------------------------------------------------------------

def ordercb(param):
    if type(param.rcb[0]) == numpy.float64:
        return param
    else:
        for i in range(len(param)):
            param.Mcb[i] = sorted(param.Mcb[i])
            param.rcb[i] = sorted(param.rcb[i])
            param.Pcb[i] = sorted(param.Pcb[i], reverse = True)
            param.Tcb[i] = sorted(param.Tcb[i], reverse = True)
            param.Ucb[i] = sorted(param.Ucb[i])
            param.Egcb[i] = sorted(param.Egcb[i], reverse = True)
            param.Etotcb[i] = sorted(param.Etotcb[i], reverse = True)
        
        return param


def aligncb(param, keyword):

    #ordercb(param)

    if keyword == 'rcb':
        val = param.rcb
    if keyword == 'Mcb':
        val = param.Mcb
    if keyword == 'Pcb':
        val = param.Pcb
    if keyword == 'Tcb':
        val = param.Tcb
    if keyword == 'Ucb':
        val = param.Ucb
    if keyword == 'Egcb':
        val = param.Egcb
    if keyword == 'Etotcb':
        val = param.Etotcb

    if type(val[0]) == numpy.float64:
        return val

    else:
        length = []
        
        for i in range(len(param)):
            length = numpy.append(length, len(val[i]))
        l = int(length.min())

        if l == 1:
            for i in range(len(param)):
                for j in range(len(val[i]) - 1):
                    val[i] = numpy.delete(val[i], 0)

        elif l == 3:

            ivec = []
            
            for i in range(len(param) - 1):

                if len(val[i]) != len(val[i + 1]):
                    ivec = numpy.append(ivec, i)
                        
            for ind in range(len(ivec)):

                if ind == 0:
                    rangei = range(int(ivec[ind]) + 1)[::-1]

                elif ind != len(ivec) - 1:
                    rangei = range(int(ivec[ind]), int(ivec[ind + 1]))

                else:
                    rangei = range(0)
                    

                for i in rangei:  #(int(ivec[0])+1)[::-1]:
                
                    if len(val[ivec[ind] + 1]) < len(val[ivec[ind]]):
                        minl = val[i + 1]
                        maxl = val[i]
                    else:
                        minl = val[i]
                        maxl = val[i + 1]                    

                    for j1 in range(3):

                        temp = []
                        for j2 in range(j1, 5):
                            temp = numpy.append(temp, abs(minl[j1] - maxl[j2]))
                        jcorr = list(temp).index(temp.min()) + j1

                        maxl[j1] = maxl[jcorr]
                        if j1 != jcorr:
                            maxl[jcorr] = 0

        for i in range(len(val)):
            #for j in range(l + 2):
            length = len(val[i])
            while length > l:
                val[i] = list(val[i])
                val[i].pop()
                length = length - 1
##                
##                try:
##                    if j >= l:
##                        val[i].pop(j)
##                except IndexError:
##                    pass
        
        return val


def aligncball(param):

    aligncb(param, 'rcb')
    aligncb(param, 'Mcb')
    aligncb(param, 'Pcb')
    aligncb(param, 'Tcb')
    aligncb(param, 'Ucb')
    aligncb(param, 'Egcb')
    aligncb(param, 'Etotcb')

    ordercb(param)

    return param

def cb_options(param):

    aligncball(param)
    
    if type(param.rcb[0]) == numpy.float64:
        return None
    else:
        length = []
        
        for i in range(len(param)):
            length = numpy.append(length, len(param.rcb[i]))
        l = int(length.min())

        if l == 1:
            return 'rcbout', 'rcbin'
        elif l == 3:
            return 'rcbout', 'rad1top', 'rcbin'
        elif l == 5:
            return 'rcbout', 'rad2top', 'rad2bottom', 'rad1top', 'rcbin'
        

def cooling_global(atmset, atmprofile, prms = prms, out = 'rcb', outrad = None, \
                   checktop = 0):
    
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
        
        if outrad == None:
                        
            Tout = atmset.Tcb
            Pout = atmset.Pcb 
            rout = atmset.rcb * Re
            Mout = atmset.Mcb * Me 
            Eout = atmset.Etotcb
            
        else:
            
            aligncball(atmset)
        
            Tout, Pout, rout, Mout, Eout = [], [], [], [], []

            if outrad == 'rcbout':
                for i in range(n):
                    Tout = numpy.append(Tout, atmset.Tcb[i][-1])
                    Pout = numpy.append(Pout, atmset.Pcb[i][-1])
                    rout = numpy.append(rout, atmset.rcb[i][-1] * Re)
                    Mout = numpy.append(Mout, atmset.Mcb[i][-1] * Me)
                    Eout = numpy.append(Eout, atmset.Etotcb[i][-1])

            elif outrad == 'rcbin':
                for i in range(n):
                    Tout = numpy.append(Tout, atmset.Tcb[i][0])
                    Pout = numpy.append(Pout, atmset.Pcb[i][0])
                    rout = numpy.append(rout, atmset.rcb[i][0] * Re)
                    Mout = numpy.append(Mout, atmset.Mcb[i][0] * Me)
                    Eout = numpy.append(Eout, atmset.Etotcb[i][0])

            elif outrad == 'rad1top':
                for i in range(n):
                    Tout = numpy.append(Tout, atmset.Tcb[i][1])
                    Pout = numpy.append(Pout, atmset.Pcb[i][1])
                    rout = numpy.append(rout, atmset.rcb[i][1] * Re)
                    Mout = numpy.append(Mout, atmset.Mcb[i][1] * Me)
                    Eout = numpy.append(Eout, atmset.Etotcb[i][1])

            elif outrad == 'rad2bottom':
                for i in range(n):
                    Tout = numpy.append(Tout, atmset.Tcb[i][2])
                    Pout = numpy.append(Pout, atmset.Pcb[i][2])
                    rout = numpy.append(rout, atmset.rcb[i][2] * Re)
                    Mout = numpy.append(Mout, atmset.Mcb[i][2] * Me)
                    Eout = numpy.append(Eout, atmset.Etotcb[i][2])

            elif outrad == 'rad2top':
                for i in range(n):
                    Tout = numpy.append(Tout, atmset.Tcb[i][3])
                    Pout = numpy.append(Pout, atmset.Pcb[i][3])
                    rout = numpy.append(rout, atmset.rcb[i][3] * Re)
                    Mout = numpy.append(Mout, atmset.Mcb[i][3] * Me)
                    Eout = numpy.append(Eout, atmset.Etotcb[i][3])                    

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

    if checktop == 0 or checktop == -1:
            
        uout = 10**interplog10u(numpy.log10(Tout), numpy.log10(Pout))
        eaccout = uout - G * Mout / rout
        deltamout = Mout[1:] - Mout[:-1]
        deltae = Eout[1:] - Eout[:-1]
        eaccav = (eaccout[1:] + eaccout[:-1]) / 2
        Lav = (L[1:] + L[:-1]) / 2
        Poutav = (Pout[1:] + Pout[:-1]) / 2

    else:
        Cv = Cvfn(prms.Y, 2./7)
        uout = Cv * Tout
        

    
    
    for i in range(len(L) - 1):
        
        if deltamout[i] > 0 and Mout[i] != 0:
            
            mnewprof = atmprofile.m[i+1] #mass profile of 'new' atmosphere (i+1)
            rnewprof = atmprofile.r[i+1]#radius profile of 'new' atmosphere(i+1)
            f = interp1d(mnewprof, rnewprof) #interpolation function
            routold = f(Mout[i]) #we find the radius in the new atmosphere (i+1)
                                 #at which the mass is equal to the mass of the
                                 #old atmosphere (i) in order to be able to
                                 #calculate dV at constant m
            deltav[i] = (4 * numpy.pi / 3) * (routold**3 - rout[i]**3)
            
        else:
            pass

        
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

    

    if type(param.rcb[0]) == numpy.float64:
        Mcb = param.Mcb * Me
        Mcbav = (Mcb[1:] + Mcb[:-1]) / 2
    else:
        param = aligncb(param)

        if len(param.rcb[0]) == 3:
            Mcbavin, Mcbavin2, Mcbavout = [], [], []
            for i in range(len(param) - 1):
                Mcbavin = numpy.append(Mcbavin, Me * (param.Mcb[i][0] + param.Mcb[i+1][0]) / 2) 
                Mcbavin2 = numpy.append(Mcbavin, Me * (param.Mcb[i][1] + param.Mcb[i+1][1]) / 2) 
                Mcbavout = numpy.append(Mcbavin, Me * (param.Mcb[i][2] + param.Mcb[i+1][2]) / 2) 

        elif len(param.rcb[0]) == 5:
            Mcbavin, Mcbavin2, Mcbavin3, Mcbavin4, Mcbavout = [], [], [], [], []
            for i in range(len(param) - 1):
                Mcbavin = numpy.append(Mcbavin, Me * (param.Mcb[i][0] + param.Mcb[i+1][0]) / 2)
                Mcbavin2 = numpy.append(Mcbavin, Me * (param.Mcb[i][1] + param.Mcb[i+1][1]) / 2) 
                Mcbavin3 = numpy.append(Mcbavin, Me * (param.Mcb[i][2] + param.Mcb[i+1][2]) / 2)
                Mcbavin4 = numpy.append(Mcbavin, Me * (param.Mcb[i][3] + param.Mcb[i+1][3]) / 2) 
                Mcbavout = numpy.append(Mcbavin, Me * (param.Mcb[i][4] + param.Mcb[i+1][4]) / 2) 
                
    Mtot = param.Mtot * Me
    MB = param.MB * Me
    
    if out == 'rcb':
        if type(param.rcb[0]) == numpy.float64:
            M = Mcb
        else:
            param = aligncb(param)
            
            if len(param.rcb[0]) == 3 or len(param.rcb[0]) == 5: 
                M = []
                for i in range(len(param)):
                    M = numpy.append(M, param.Mcb[i][-1] * Me)

    elif out == 'RHill':
        M = Mtot
    elif out == 'RB':
        M = MB
        
    Mav = (M[1:] + M[:-1]) / 2
    Ldt = 0 * numpy.ndarray(shape = (n - 1), dtype = float)
    Ldtr = 0 * numpy.ndarray(shape = (n - 1), dtype = float)

    #if type(param.rcb[0]) == numpy.float64:
    
    for i in range(n - 1):

        if M[i + 1] > M[i]:

            if onlyrad == 0:
                m = numpy.linspace(prms.Mco, M[i], npoints)
            elif onlyrad != 0 and out == 'RHill':
                if type(param.rcb[0]) == numpy.float64:
                    m = numpy.linspace(Mcbav[i], Mtot[i], npoints)
                else:
                    param = aligncb(param)
                    if len(param.rcb[0]) == 3:
                        m = numpy.append(numpy.linspace(Mcbavin[i], Mcbavin2[i], npoints), \
                                         numpy.linspace(Mcbavout[i], Mtot[i], npoints))
                    elif len(param.rcb[0]) == 5:
                        m = numpy.append( \
                            numpy.append(numpy.linspace(Mcbavin[i], Mcbavin2[i], npoints), \
                                         numpy.linspace(Mcbavin3[i], Mcbavin4[i], npoints)), \
                                         numpy.linspace(Mcbavout[i], Mtot[i], npoints))
                    
                        
            elif onlyrad != 0 and out == 'RB':
                if type(param.rcb[0]) == numpy.float64:
                    m = numpy.linspace(Mcbav[i], MB[i], npoints)
                else:
                    param = aligncb(param)
                    if len(param.rcb[0]) == 3:
                        m = numpy.append(numpy.linspace(Mcbavin[i], Mcbavin2[i], npoints), \
                                         numpy.linspace(Mcbavout[i], MB[i], npoints))
                    elif len(param.rcb[0]) == 5:
                        m = numpy.append( \
                            numpy.append(numpy.linspace(Mcbavin[i], Mcbavin2[i], npoints), \
                                         numpy.linspace(Mcbavin3[i], Mcbavin4[i], npoints)), \
                                         numpy.linspace(Mcbavout[i], MB[i], npoints))
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

            dS = 10**interplog10S(numpy.log10(T2int), numpy.log10(P2int)) - \
                 10**interplog10S(numpy.log10(T1int), numpy.log10(P1int))

            Ldt[i] = - sum(Tav * dS * dm)
            
        else:
            i = i + 1

    return Ldt

##    else:
##
##        for i in range(n - 1):
##
##            if M[i + 1] > M[i]:
##
##                if onlyrad == 0:
##                    m = numpy.linspace(prms.Mco, M[i], npoints)
##                elif onlyrad != 0 and out == 'RHill':
##                    m = numpy.linspace(Mcbav[i][0], Mtot[i], npoints)
##                    #mr = numpy.linspace(Mcbav[i][-1], Mtot[i], npoints)
##                elif onlyrad != 0 and out == 'RB':
##                    m = numpy.linspace(Mcbav[i][0], MB[i], npoints)
##                    #mr = numpy.linspace(Mcbav[i][-1], MB[i], npoints)
##                else:
##                    print "Wrong choice of boundary."
##                    sys.exit()
##
##                mmid = (m[:-1] + m[1:]) / 2
##                dm = m[1] - m[0]
##                
##                m1 = prof.m[i]
##                T1 = prof.t[i]
##                P1 = prof.P[i]
##
##                m2 = prof.m[i + 1]
##                T2 = prof.t[i + 1]
##                P2 = prof.P[i + 1]
##
##                fP1 = interp1d(m1, P1)
##                fT1 = interp1d(m1, T1)
##                
##                fP2 = interp1d(m2, P2)
##                fT2 = interp1d(m2, T2)
##
##                success = 0
##            
##                while success != 1:
##                    try:
##                        P1int = fP1(mmid)
##                        P2int = fP2(mmid)
##                        T1int = fT1(mmid)
##                        T2int = fT2(mmid)
##
##                        success = 1
##
##                    except ValueError:
##                        mmid = mmid[1:]
##
##                Tav = (T1int + T2int) / 2
##
##                dS = 10**interplog10S(numpy.log10(T2int), numpy.log10(P2int)) - \
##                     10**interplog10S(numpy.log10(T1int), numpy.log10(P1int))
##
##                Ldt[i] = - sum(Tav * dS * dm)
##
####                mmidr = (mr[:-1] + mr[1:]) / 2
####                dmr = mr[1] - mr[0]
####
####                success = 0
####            
####                while success != 1:
####                    try:
####                        P1intr = fP1(mmidr)
####                        P2intr = fP2(mmidr)
####                        T1intr = fT1(mmidr)
####                        T2intr = fT2(mmidr)
####
####                        success = 1
####
####                    except ValueError:
####                        mmidr = mmidr[1:]
####
####                Tavr = (T1intr + T2intr) / 2
####
####                dSr = 10**interplog10S(numpy.log10(T2intr), numpy.log10(P2intr)) - \
####                     10**interplog10S(numpy.log10(T1intr), numpy.log10(P1intr))
####
####                Ldtr[i] = - sum(Tavr * dSr * dmr)
##                
##            else:
##                i = i + 1
##
##        return Ldt #, Ldtr, Ldt + Ldtr

        
    

def critical(param, prof, prms = prms, out = 'rcb', outrad = None, checktop = 0):
    
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

    param = aligncball(param)
    dt = cooling_global(param, prof, prms, out = out, outrad = outrad, checktop = checktop)[0]
    Lav = (param.L[:-1] + param.L[1:]) / 2

    #if type(param.rcb[0]) == numpy.float64:
    
    for i in range(len(dt)):
        if dt[i] < 0:
            break
    for j in range(len(param.MB)):
        if param.MB[j] > 2 * prms.Mco / Me:
            break
        
    k = numpy.array([i, j]).min()
    return param[:k], prof[:k], dt[:k]


def Lrad(param, prof, n):

    #n = len(param)
    
    fP = interp1d(prof[n].r, prof[n].P)
    frho = interp1d(prof[n].r, prof[n].rho)
    fm = interp1d(prof[n].r, prof[n].m)
    fT = interp1d(prof[n].r, prof[n].t)
    fdelad = interp1d(prof[n].r, prof[n].delad)

    P = fP(param.rcb[n] * Re)
    #Pf = fP(param.rcb * Re)

    rho = frho(param.rcb[n] * Re)
    #rhof = frho(rcball[n,1] * Re)

    m = fm(param.rcb[n] * Re)
    #mf = fm(rcball[n,1] * Re)

    T = fT(param.rcb[n] * Re)
    #T = fT(param.rcb * Re)

    delad = fdelad(param.rcb[n] * Re)
    #deladf = fdelad(rcball[n,1] * Re)

    k = interp_opacity(T)
    #kf = kdustall(Tf, rhof)

    L = 64 * numpy.pi * G * m * sigma * T**4 * delad / (3 * k * P)
    #Lf = 64 * pi * G * mf * sigma * Tf**4 * deladf / (3 * kf * Pf)
    
    return L

##def aligncb2(param):
##    
##    param = ordercb(param)
##
##    if type(param.rcb[0]) == numpy.float64:
##        return param
##    else:
##        length = []
##        for i in range(len(param)):
##            length = numpy.append(length, len(param.rcb[i]))
##        l = int(length.max())
##
##        for i in range(len(param)):
##            if len(param.rcb[i]) == l:
##                pass
##            else:
##                dl = l - len(param.rcb[i])
##                for j in range(dl):
##                    param.rcb[i] = numpy.append(0, param.rcb[i])
##        return param
##
##
##def aligncb_new(param):
##    
##    param = ordercb(param)
##
##    if type(param.rcb[0]) == numpy.float64:
##        return param
##
##    else:
##        ivec = []
##
##        for i in range(len(param) - 1):
##            if len(param.rcb[i]) != len(param.rcb[i + 1]):
##                ivec = numpy.append(ivec, i)
##
##        for i in range(len(ivec)):
##            l = len(param.rcb[ivec[i] + 1]) - len(param.rcb[ivec[i]])
##            
##            if l > 0:
##                minl = param.rcb[ivec[i]]
##                maxl = param.rcb[ivec[i] + 1]
##            else:
##                minl = param.rcb[ivec[i] + 1]
##                maxl = param.rcb[ivec[i]]
##
##            l = numpy.abs(l)
##
##            for j in range(l):
##                minl = numpy.append(0, minl)
##
##            #diffmatrix = 0 * numpy.ndarray(shape = (len(minl), len(maxl)), \
##            #                               dtype = float)
##            for imin in range(len(minl) - 1):
##                diff = []
##                for imax in range(len(maxl) - 1):
##                    diff = numpy.append(diff, numpy.abs(minl[imin] - maxl[imax]))

                        
##                                param.rcb[k] = minl



##
##
##                            
##                            lp = len(param.rcb[i + 1]) - len(param.rcb[i])
##                
##                    if lp < 0:
##                        for k in range(i + 2)[::-1]:
##                        for j in range(-l):
##                        
##                            param.rcb[k] = numpy.delete(param.rcb[k], 0)
                    

##def aligncb(param):
##    
##    ordercb(param)
##
##    if type(param.rcb[0]) == numpy.float64:
##        return param
##
##    else:
##        length = []
##        
##        for i in range(len(param)):
##            length = numpy.append(length, len(param.rcb[i]))
##        l = int(length.min())
##
##        if l == 1:
##            for i in range(len(param)):
##                for j in range(len(param.rcb[i]) - 1):
##                    param.rcb[i] = numpy.delete(param.rcb[i], 0)
##
##        elif l == 3:
##            ivec = []
##
##            for i in range(len(param) - 1):
##                if len(param.rcb[i]) != len(param.rcb[i + 1]):
##                    ivec = numpy.append(ivec, i)
##
##            for i in range(len(ivec)):
##                
##                if len(param.rcb[ivec[i] + 1]) > len(param.rcb[ivec[i]]):
##                    #minl = param.rcb[ivec[i]]
##                    #maxl = param.rcb[ivec[i] + 1]
##
##                    if len(ivec) == 1:
##                        
##                        for k in range(int(ivec[i]))[::-1]:
##                                      
##                            minl = param.rcb[k]
##                            maxl = param.rcb[k + 1]
##
##                            for j in range(2):
##                                minl = numpy.append(0, minl)
##
##                            for j1 in range(2, 4):
##                                temp = []
##                                for j2 in range(4):
##                                    temp = numpy.append(temp, abs(minl[j1] - maxl[j2]))
##                                jcorr = list(temp).index(temp.min())
##
##                                minl[jcorr] = minl[j1]
##                                if j1 != jcorr:
##                                    minl[j1] = 0
##                                param.rcb[k] = minl
##
##                    else:
##
##                        for k in range(int(ivec[i]), int(ivec[i + 1])):
##                                      
##                            minl = param.rcb[k + 1]
##                            maxl = param.rcb[k]
##
##                            for j in range(2):
##                                minl = numpy.append(0, minl)
##
##                            for j1 in range(2, 4):
##                                temp = []
##                                for j2 in range(4):
##                                    temp = numpy.append(temp, abs(minl[j1] - maxl[j2]))
##                                jcorr = list(temp).index(temp.min())
##
##                                minl[jcorr] = minl[j1]
##                                if j1 != jcorr:
##                                    minl[j1] = 0
##                                param.rcb[k + 1] = minl
##
##
##                else:
##                    if len(ivec) == 1:
##                        
##                        for k in range(int(ivec[i]), len(param) - 1):
##                                      
##                            minl = param.rcb[k + 1]
##                            maxl = param.rcb[k]
##
##                            for j in range(2):
##                                minl = numpy.append(0, minl)
##
##                            for j1 in range(2, 4):
##                                temp = []
##                                for j2 in range(4):
##                                    temp = numpy.append(temp, abs(minl[j1] - maxl[j2]))
##                                jcorr = list(temp).index(temp.min())
##
##                                minl[jcorr] = minl[j1]
##                                if j1 != jcorr:
##                                    minl[j1] = 0
##                                param.rcb[k + 1] = minl
##
##                    else:
##
##                        for k in range(int(ivec[i]), int(ivec[i + 1])):
##                                      
##                            minl = param.rcb[k + 1]
##                            maxl = param.rcb[k]
##
##                            for j in range(2):
##                                minl = numpy.append(0, minl)
##
##                            for j1 in range(2, 4):
##                                temp = []
##                                for j2 in range(4):
##                                    temp = numpy.append(temp, abs(minl[j1] - maxl[j2]))
##                                jcorr = list(temp).index(temp.min())
##
##                                minl[jcorr] = minl[j1]
##                                if j1 != jcorr:
##                                    minl[j1] = 0
##                                param.rcb[k + 1] = minl
##
##        return param


##def insert_zeros(param):
##
##    aligncb(param)
##
##    if type(param.rcb[0]) == numpy.float64:
##        return param
##
##    else:
##        
##        for i in range(len(param)):
##            for j in range(len(param.rcb[i])):
##                if param.rcb[i][j] == 0:
##                    for k in range(len(param)):
##                        param.rcb[k][j] = 0 #list(param.rcb[k])
##                        #param.rcb[k].pop(j)
##        return param

##        jvec = []                
##        for j in range(len(param.rcb[5])):
##            if param.rcb[5][j] == 0:
##                jvec = numpy.append(jvec, j)
##                
##        for i in range(len(param)):
##            param.rcb[i] = list(param.rcb[i])
##            for k in range(len(jvec)):
##                temp = param.rcb[i]
##                param.rcb[i].pop(int(jvec[k]))
##        
##        return param




    
                    
        

##def remove_zeros(param):
##
##    param = aligncb2(param)
##
##    if type(param.rcb[0]) == numpy.float64:
##        return param
##    
##    else:
##        n = len(param)
##        ncol = len(param.rcb[0])
##        
##        for j in range(ncol - 1):
##            for i in range(n):
##                if param.rcb[i][j] == 0:
##                    if param.rcb
##                    jplus, jminus = j, j
##                    while param.rcb[i][jplus] == 0 and param.rcb[i][jminus] == 0:
##                        jplus = jplus + 1
##                        jminus = jminus - 1
##                    if param.rcb[i][jplus] != 0:
##                        
##                    
##                    
##                    #for k in range(i + 1, n):
##                    #    if param.rcb[k][j] != 0:
##                    #    break
##                #if k != n-1:
##                    
##                    
##
##        
####        length = []
####        for i in range(len(param)):
####            length = numpy.append(length, len(param.rcb[i]))
####        l = int(length.max())
####
##        ivec = []
##
##        for i in range(len(param) - 1):
##            if len(param.rcb[i]) != len(param.rcb[i + 1]):
##                ivec = numpy.append(ivec, i)
##                
##        for i in range(len(ivec)):
##            if len(param.rcb[ivec[i]]) <= len(param.rcb[ivec[i] + 1]):
####                minl = param.rcb[ivec[i]]
####                maxl = param.rcb[ivec[i] + 1]
####            else:
####                minl = param.rcb[ivec[i] + 1]
####                maxl = param.rcb[ivec[i]]
####                
####            for j1 in range(len(minl)):
####                    temp = []
####                    for j2 in range(len(maxl)):
####                        temp = numpy.append(temp, numpy.abs( \
####                            minl[j1] - maxl[j2]))
####                    correct = temp.min()   
####                    for j2 in range(len(maxl)):
####                        if numpy.abs(minl[j1] - maxl[j2]) != correct:
####                            maxl[j2] = 0
##        return param
##                                        
##
##        
####        for j in range(len(param.rcb[i])):
####            if param.rcb[i][j] == 0:
####                pass
####            else:
####                for k in range(i)[::-1]:
####                    if param.rcb[k]
##            
                
            



##def aligncb(param):
##    param = ordercb(param)
##
##    if type(param.rcb[0]) == numpy.float64:
##        return param
##    
##    else:
##        for i in range(len(param) - 1):
##            if len(param.rcb[i]) == len(param.rcb[i + 1]):
##                pass
##            else:
##                l = len(param.rcb[i + 1]) - len(param.rcb[i])
##                
##                if l < 0:
##                    for k in range(i + 1):
##                        for j in range(-l):
##                        
##                            param.rcb[k] = numpy.delete(param.rcb[k], 0)
##                            param.Pcb[k] = numpy.delete(param.Pcb[k], 0)
##                            param.Tcb[k] = numpy.delete(param.Tcb[k], 0)
##                            param.Mcb[k] = numpy.delete(param.Mcb[k], 0)
##                            param.Ucb[k] = numpy.delete(param.Ucb[k], 0)
##                            param.Egcb[k] = numpy.delete(param.Egcb[k], 0)
##                            param.Etotcb[k] = numpy.delete(param.Etotcb[k], 0)
##                        
##                else:
##                    
##                    for j in range(l):
##                        
##                        param.rcb[i + 1] = numpy.delete(param.rcb[i + 1], 0)
##                        param.Pcb[i + 1] = numpy.delete(param.Pcb[i + 1], 0)
##                        param.Tcb[i + 1] = numpy.delete(param.Tcb[i + 1], 0)
##                        param.Mcb[i + 1] = numpy.delete(param.Mcb[i + 1], 0)
##                        param.Ucb[i + 1] = numpy.delete(param.Ucb[i + 1], 0)
##                        param.Egcb[i + 1] = numpy.delete(param.Egcb[i + 1], 0)
##                        param.Etotcb[i + 1] = numpy.delete(param.Etotcb[i + 1], 0)
##                        
##        return param

