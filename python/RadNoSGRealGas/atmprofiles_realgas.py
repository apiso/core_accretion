"""

This module calculates the minimum entropy of an atmosphere set for given
initial conditions, and generates atmosphere profiles & plots.

"""


from utils.constants import G, kb, mp, Rb, Me, Re, RH, RHe, sigma
import numpy
import scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate
from utils.table_size import Nt, Np, Ncol
from utils.parameters import rc, Mc, Y, beta, R, gamma, delad, Cv, a
from utils.disk_model import Td, Pd, rhod, kappa
import utils.cb_corrections as cb
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from shooting_realgas import shoot, Tcore
from utils.interpolation_functions_logp import interpolationdelad


#---------------------------------------------------------------------------


def atmhelp(x, S, T1, T2, n):
    """
    Returns atmosphere profile for given S, Td, Pd, Mc, rc:
    
    Returns an array of the follwing:
        T
        r(T) (array)
        P(T) (array)
        rho(T) (array)
        u(T) (array)
        m(T) (array)
        M_ad = mass of convective region
        M_iso = mass of radiative region
        U = internal energy of convective region
        E = gravitational energy of the convective region
        U + E = total energy of convective region
        U_iso = internal energy of radiative region
        Eg_iso = gravitational energy of radiative region
        E_iso = total energy of radiative region
        fact = P/(rho * delad)
        RB = Bondi radius
        RHill = Hill radius
        L = luminosity in the convective region
    """    
    
    t = x[0]
    r = x[1]
    m = x[2]
    P = x[3]
    rho = x[4]
    U = x[5] 
    E = x[6]
    vir = x[7]
    RB = x[-2]
    RHill = x[-1]

    rb = r[-1] #radius of convective zone
    Mcb = m[-1] #mass of radiative zone
    
    def dmiso(r):
        return 4 * numpy.pi * rhod * numpy.exp(RB * (1. / r - 1. / RHill)) * r**2 #dmiso = (4*pi*r**2)*(rhod*exp(RB/r))*dr
    Miso = integrate.quad ( dmiso, rb, RHill )[0] #integrate dmiso

    Rout = RHill #(R / Cv) * RB #zero energy point
    
    def duiso(r):
        return 4 * numpy.pi * (Cv * Td) * rhod * numpy.exp(RB * (1. / r - 1. / RHill)) * r**2
        #duiso = Cv*Td*dm = Cv*Td*(4*pi*r**2*rhod*exp(Rb/r))

    Uiso = integrate.quad(duiso, rb, Rout)[0] #integrate dUiso 


    def degiso(r):
        return - G * m[-1] * 4 * numpy.pi * r * rhod * \
               numpy.exp(RB * (1. / r - 1. / RHill))
    Egiso = integrate.quad(degiso, rb, Rout)[0]
    Eiso = Uiso + Egiso

    L = 64 * numpy.pi * G * m[-1] * sigma * t[-1]**4 * interpolationdelad(numpy.log10(P[-1]), S) / (3 * kappa(t[-1]) * P[-1])

           
    return numpy.array ([t, r, m, P, rho, Mcb / Me, Miso / Me, rb / Re, RB / Re, RHill / Re, \
                         P[0], P[-1], t[0], t[-1], E, U, E + U, Egiso, Uiso, Eiso, vir, L])



def atm(S, T1, T2, n):
    x = shoot(S, T1, T2, n, 10**(-4))
    x1 = x[0]
    x2 = x[1]
    atm1 = atmhelp(x1, S, T1, T2, n)
    atm2 = atmhelp(x2, S, T1, T2, n)
    return atm1, atm2



def atmparam(S, T1, T2, n):

    """
    Returns Tc, Pc, Pb, rb, RB, M_convective, M_radiative, U, E_grav, E_tot, E_iso, virial check
    for a model atmosphere
    """
    
    x1 = atm(S, T1, T2, n)[0]
    x2 = atm(S, T1, T2, n)[1]
    param1 = numpy.array([x1[0][0], x1[3][0], x1[3][-1], x1[1][-1]/Re, x1[-1]/Re, x1[5], \
                         x1[6], x1[7], x1[8], x1[9], x1[12], x1[13]])
    param2 = numpy.array([x2[0][0], x2[3][0], x2[3][-1], x2[1][-1]/Re, x2[-1]/Re, x2[5], \
                         x2[6], x2[7], x2[8], x2[9], x2[12], x2[13]])
    return param




def atmtot(S1, S2, npoints, n, T1, T2):

    """
    Returns a set of atmosphere profiles and parameters for a series on entropies.

    Input:
        - S1 = minimum entropy
        - S2 = maximum entropy
        - npoints = number of entropy points (integer)
        - n = number of grid points
        - T1 = first temperature guess 
        - T2 = second temperature guess

    Output:
        - two 3D arrays (profile1 for low mass case, profile2 for high mass case) of atmosphere
        profiles: for each entropy, it returns t, r(t), m(t), p(t), rho(t)
        - two 2D arrays (param1 for low mass case, param2 for high mass case) of atmosphere
        parameters: for each entropy, ir returns Tc, Pc, Pb, rb, RB, RHill,
        M_convective, M_radiative, U, E_grav, E_tot, E_iso, virial check, L
    """

    S = numpy.linspace(S1, S2, npoints)
    profile1 = 0*numpy.ndarray(shape = (len(S), n, 5), dtype = float ) #initializing empty arrays
    profile2 = 0*numpy.ndarray(shape = (len(S), n, 5), dtype = float ) #for atm profile
    param1 = 0*numpy.ndarray(shape = (len(S), 17), dtype = float ) #initializing empty arrays 
    param2 = 0*numpy.ndarray(shape = (len(S), 17), dtype = float ) #for atm parameters
        
        
    for i in range(len(S)):
        x = atm(S[i], T1, T2, n)
        x1 = x[0]
        x2 = x[1]
        for j in range(len(x1[0])):
            profile1[i,j,0] = x1[0][j]
            profile1[i,j,1] = x1[1][j]
            profile1[i,j,2] = x1[2][j]
            profile1[i,j,3] = x1[3][j]
            profile1[i,j,4] = x1[4][j]
            
        param1[i] = numpy.array([x1[5], x1[6], x1[7], x1[8], x1[9], x1[10], x1[11], x1[12], x1[13], \
                                x1[14], x1[15], x1[16], x1[17], x1[18], x1[19], x1[20], x1[21]])
        
        for j in range(len(x2[0])):
            profile2[i,j,0] = x2[0][j]
            profile2[i,j,1] = x2[1][j]
            profile2[i,j,2] = x2[2][j]
            profile2[i,j,3] = x2[3][j]
            profile2[i,j,4] = x2[4][j]
            
        param2[i] = numpy.array([x2[5], x2[6], x2[7], x2[8], x2[9], x2[10], x2[11], x2[12], x2[13], \
                                x2[14], x2[15], x2[16], x2[17], x2[18], x2[19], x2[20], x2[21]])
    return profile1, param1, profile2, param2



#####################################################################################################
"""
Various plots for model atmospheres
"""



def ptplot(atm):
    
    """P-T profile"""
    
    rnb = numpy.linspace(0, 255, numpy.shape(atm)[0])
    for i in range(numpy.shape(atm)[0]):
        plt.loglog(atm[i,:,3], atm[i,:,0], linewidth=2.0, \
                   c=cm.rainbow(numpy.int(rnb[i])))
    plt.ylabel('T (K)')
    plt.xlabel(r'P (dyn cm$^{-2}$)')
    plt.title(r'a=%s AU, $M_c=$%s$M_{\oplus}$, $Y=$%s ' %(str(a), str(Mc/Me), str(Y)))

def rtplot(atm):

    """r-T profile"""
    
    rnb = numpy.linspace(0, 256, numpy.shape(atm)[0])
    for i in range(numpy.shape(atm)[0]):
        plt.loglog(atm[i,:,0], atm[i,:,1]/rc, linewidth=2.0, \
                   c=cm.rainbow(numpy.int(rnb[i])))
    plt.xlabel('log(T)')
    plt.ylabel(r'r ($r_c$)')
    plt.title(r'a=%s AU, $M_c=$%s$M_{\oplus}$, $Y=$%s ' %(str(a), str(Mc/Me), str(Y)))
    return

def mtplot(atm):

    """M-T profile"""
    
    rnb = numpy.linspace(0, 255, numpy.shape(atm)[0])
    for i in range(numpy.shape(atm)[0]):
        plt.loglog(atm[i,:,2]/Me, atm[i,:,0], linewidth=2.0, \
                   c=cm.rainbow(numpy.int(rnb[i])))
    plt.ylabel('T (K)')
    plt.xlabel(r'$M_{\mathrm{atm}}$ ($M_{\oplus}$)')
    plt.title(r'a=%s AU, $M_c=$%s$M_{\oplus}$, $Y=$%s ' %(str(a), str(Mc/Me), str(Y)))

def mrplot(atm):

    """M-r profile"""
    
    rnb = numpy.linspace(0, 255, numpy.shape(atm)[0])
    for i in range(numpy.shape(atm)[0]):
        plt.loglog(atm[i,1:,2]/Me, atm[i,1:,1]/Re, linewidth=2.0, \
                   c=cm.rainbow(numpy.int(rnb[i])))
    plt.ylabel(r'r ($r_{\oplus}$)')
    plt.xlabel(r'$M_{\mathrm{atm}}$ ($M_{\oplus}$)')
    plt.title(r'a=%s AU, $M_c=$%s$M_{\oplus}$, $Y=$%s ' %(str(a), str(Mc/Me), str(Y)))

def mrhoplot(atm):

    """rho-M profile"""
    
    rnb = numpy.linspace(0, 255, numpy.shape(atm)[0])
    for i in range(numpy.shape(atm)[0]):
        plt.loglog(atm[i,:,2]/Me, atm[i,:,4], linewidth=2.0, \
                   c=cm.rainbow(numpy.int(rnb[i])))
    plt.ylabel(r'$\rho$ (g cm$^{-3}$)')
    plt.xlabel(r'$M_{\mathrm{atm}}$ ($M_{\oplus}$)')
    plt.title(r'a=%s AU, $M_c=$%s$M_{\oplus}$, $Y=$%s ' %(str(a), str(Mc/Me), str(Y)))

def lmplot(atmset):
    """Plots luminosity vs. mass of the convective zone for a given atmosphere set"""
    plt.semilogy(atmset[:, 0], atmset[:, -1])
    plt.xlabel(r'$M_{\mathrm{atm}}$($M_{\oplus}$)')
    plt.ylabel(r'L (erg s$^{-1}$)')
    plt.title(r'a=%s AU, $M_c=$%s$M_{\oplus}$, $Y=$%s ' %(str(a), str(Mc/Me), str(Y)))

def smplot(atmset, S):
    """Plots mass of convective zone vs. entropy for a given atmosphere set"""
    plt.semilogx(atmset[:, 0], S)
    plt.ylabel('log(S)')
    plt.xlabel(r'$M_{\mathrm{atm}}$ ($M_{\oplus}$)')
    plt.title(r'a=%s AU, $M_c=$%s$M_{\oplus}$, $Y=$%s ' %(str(a), str(Mc/Me), str(Y)))

#################################################################################

##
##def mins(n, T1, T2, tol, tolt, S, step):
##
##    """
##    Finds the minimum entropy of a given atmosphere set.
##    """
##    
##    #S = interpolations(Td, Pd, T, EOS) - 0.001
##    Tc1 = T1
##    Tc2 = T2
##
##    sol = Tcore(S, min(Tc1, Tc2), T2, n, tol)
##    Tc1 = sol[0]
##    Tc2 = sol[1]
##    
##    while (sol != 1):
##        Smax = S
##        print "S is %s" %S
##        S = S - step
##        sol = Tcore(S, min(Tc1, Tc2), T2, n, tol)
##        if sol != 1:
##            Tc1 = sol[0]
##            Tc2 = sol[1]
##    
##    Slow = S
##    S = (S + Smax) / 2
##
##    while (abs(Tc1 - Tc2) > tolt):
##        sol = Tcore(S, min(Tc1, Tc2), T2, n, tol)
##        if (sol == 1):
##            Slow = S
##            print "S is %s" %S
##            S = (S + Smax) / 2
##            #sol = Tcore(S, min(Tc1, Tc2), T2, n, tol)
##        else:
##            #if sol[0] < Tc1:
##            #    T2 = sol[0]
##            #if sol[1] < Tc2:
##            #    T2 = sol[1]
##            Tc1 = sol[0]
##            Tc2 = sol[1]
##            Smax = S
##            print "S is %s" %S
##            S = (S + Slow) / 2
##            #sol = Tcore(S, min(Tc1, Tc2), T2, n, tol)
##
##    print "S is %s" %S
##    print "Tc low is %s, Tc high is %s" %(Tc1, Tc2)
##    return Smax, Tc1, Tc2










