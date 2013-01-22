"""

Shooting code for polytrope. For a given entropy, 'Tcore' shoots for the core temperature for which mass and radius
match at the radiative convective boundary. Once the solution is found, 'shoot' generates the atmosphere profile.

AY: a bit more specific...
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


#-----------------------------------------------------------------------------

def Tcore(log10Sc, T1, T2, n, tol ):

    """
    Implements the shooting method to solve for a model adiabatic atmosphere. It shoots for the core
    temperature for an atmosphere with a given entropy. The atmosphere and 
    disk parameters (Mc, rc, Td, Pd) are taken from the parameters file.

    Input:
        log10Sc = log10 of desired entropy
        T1 = first try for core temperature
        T2 = second try for core temperature
        n = number of temperature grid points
        tol = allowable tolerance at the boundary

    Output:
        Tc1 = core temperature for the first solution
        Tc2 = core temperature for the second solution

    """

    Tb = cb.cbcorr(log10Sc)[0]  #temperature at top of convective zone
    Pb = cb.cbcorr(log10Sc)[1]  #pressure at top of convective zone
    rhob = Pb / (R * Tb) #density at top of convective zone
    k = Pb / rhob**gamma #Pb = k * rhob**gamma
    deladcb = delad #cb.delad(Tb, log10Sc) #delad at top of convective zone
    delo = deladcb * ((Tb/Td)**(4-beta) / (Pb/Pd)) #delzero
    delinf = 1. /(4 - beta)

    def g(p):
        """
        Hydrostatic balance in the radiative zone; it reduces to 1/rb - 1/RHill = (R*Td/(G*M))int_Pb^Pd{g(p)} 
        """
        return - (1./p) * (1 + (delo / delinf) * (p / Pd - 1))**(1. /(4 - beta))
    
    
    pint = integrate.quad(g, Pb, Pd)[0] #int_Pb^Pd{g(p)} 

    if (pint > 0):

        def f(x, logp):
            """
            structure eqns. "x" = [ r , m ]
            dr/dlogp = - r**2/(G * m) * k**(1/gamma) * p**(1-1/gamma),
            dm/dlogp = -4 * pi/(G * m) * p * r**4
            """
            p = numpy.exp(logp)
            return numpy.array([ - x[0]**2 /(G * x[1]) * k**(1./gamma) * p**(1 - 1/gamma), \
                                 - 4 * numpy.pi * p * x[0]**4 / ( G * x[1]) ]) 
               
            

        def delta(temp):
            """
            Returns the relative matching error for the radius of the convective zone
            """

            pinit = (R**gamma / k)**(1./(gamma - 1)) * temp**(gamma/(gamma - 1))
                                #finds the corresponding core pressure for the guessed temperature temp
            logp1 = numpy.linspace(numpy.log(pinit), numpy.log(Pb), n) #pressure grid, linear in natural log p
            #AY: why is grid needed? To find matching error just go from one end to other.
                ##AP: I am shooting for a core temperature, but the grid is in log pressure, so in order to know Pc = Ptry
                    ## I need to calculate it as a function of temp = Ttry
            
            y = odeint(f, [rc, Mc], [numpy.log(pinit), numpy.log(Pb)]) #integration of the structure equations
     
            rcb1 = y[:,0][-1] #radius of convective zone
            Mcb1 = y[:,1][-1] #mass of convective zone
            RHill1 = RHill(Mcb1, a) #Hill radius
            rcb = (1. / RHill1 + (R * Td / (G * Mcb1)) * pint)**(-1)
            deltab1 = 4 / numpy.pi * numpy.arctan(rcb1 / rcb) - 1
            if math.isnan(deltab1):
                deltab1 = 1.
            return - deltab1


        Topt, deltamax, ierr, numfunc = fminbound(delta, T1, T2, xtol = tol, full_output = 1)

        Tlow, errTlow = brentq(delta, T1, Topt, xtol = tol, full_output = 1)
        Thigh, errThigh = brentq(delta, Topt, T2, xtol = tol, full_output = 1)

        return Tlow, Thigh, delta(Tlow), delta(Thigh)
             
    else:
        print "The radius of the convective boundary is larger than the Hill radius. This is not a \
        physical solution, choose a lower entropy."
    
#-----------------------------------------------------------------------------

def shoot(log10Sc, T1, T2, n, tol):
    
    """
    Input:
        log10Sc = log10 of desired entropy
        T1 = first try for core temperature
        T2 = second try for core temperature
        n = number of temperature grid points
        tol = allowable tolerance at the boundary

    Output:
        t = array of temperatures (from Tc to Tb)
        r = array of r(t)
        m = array of m(t)
        p = array of p(t)
        rho = array of rho(t)
        atmosphere parameters (Mcb, Mrad, rcb, RB, RHill, Pc, Pcb, Tc, Tcb, Eg, U, Etot, L, vircheck)
    """
    
    Tc = Tcore(log10Sc, T1, T2, n, tol)
    Tc1 = Tc[0]
    Tc2 = Tc[1]
    err1 = Tc[2]
    err2 = Tc[3]
    Pb = cb.cbcorr(log10Sc)[1]
    Tb = cb.cbcorr(log10Sc)[0]
    rhob = Pb / (R * Tb) 
    k = Pb / rhob**gamma 
    
    def ivp(temp):
        pinit = (R**gamma / k)**(1./(gamma - 1)) * temp**(gamma/(gamma - 1))
        logp1 = numpy.linspace(numpy.log(pinit), numpy.log(Pb), n)
        p1 = numpy.exp(logp1)

        def f(x, logp):
            """
            Structure equations with virial integrals
            Returns dy/dlogp where y = (r, M, Eg, Iu), with Iu the 3p/rho dm integral in the virial theorem

            AY: I prefer a single function of structure equations with an option of also returning virial
            integrals.  Doing it that way would admittedly cause you to be more explicit about arguments passed,
            which could be a good thing.
            """
            p = numpy.exp(logp)
            
            #AY: there's no need to return numpy arrays.  A simple list will do.  Why be more complicated that necessary?
            return numpy.array([ - x[0]**2 /(G * x[1]) * k**(1./gamma) * p**(1 - 1/gamma), \
                                     - 4 * numpy.pi * p * x[0]**4 / ( G * x[1]), \
                                 4 * numpy.pi * p * x[0]**3, \
                                - 12 * numpy.pi / (G * x[1]) * x[0]**4 * k**(1./gamma) * p**(2 - 1./gamma)])

        E0 = G * Mc**2 / rc
        y = odeint(f, [rc, Mc, - E0, E0], logp1, mxstep = 50000)

        RHill1 = RHill(y[:, 1][-1], a)
        RB1 = G * (y[:, 1][-1]) / (R * Td)
        rb = y[:,0][-1]
        rho = (p1 / k)**(1./gamma)
        t = p1 / (R * rho)
        Pc = numpy.exp(logp1)[0]
        rhs = y[:,3][-1] + 4 * numpy.pi * (rc**3 * Pc - rb**3 * Pb)
        xi = 3 * (gamma - 1)
        Eg = E0 + y[:, 2][-1]
        U = (y[:, 3][-1] - E0) / xi

        r = y[:,0]
        m = y[:,1]
        P = p1
        vircheck = (rhs + y[:, 2][-1]) / abs(y[:, 2][-1])

        L = 64 * numpy.pi * G * m[-1] * sigma * t[-1]**4 * delad / (3 * kappa(t[-1]) * P[-1])


        def dmiso(r):
            return 4 * numpy.pi * rhod * numpy.exp(RB1 * (1. / r - 1. / RHill1)) * r**2
                #dmiso is computed as though radiative zone is isothermal, needs updating
        
        Miso = integrate.quad ( dmiso, rb, RHill1 )[0] / Me #integrate dmiso

        Mcb = m[-1] / Me
        rcb = r[-1]
        Pc = P[0]
        Pcb = P[-1]
        tc = t[0]

        return t, r, m, P, rho, Mcb, Miso, rcb / Re, RB1 / Re, RHill1 / Re, Pc, Pcb, tc, Tb, Eg, U, Eg + U, L, vircheck

    return (ivp(Tc1), err1), (ivp(Tc2), err2)


def minent(log10S1, log10S2, T1, T2, n, tol, tols):
    
    """
    Function to find the minimum entropy by bracketing. zbrac wasn't immediately obvious how to use,
    but this is relatively fast for good initial entropy brackets.
    """

    def helper(log10S):
        try:
            sol = Tcore(log10S, T1, T2, n, tol)
            return sol
        except ValueError:
            return (0, 10)

    sol2 = Tcore(log10S2, T1, T2, n, tol)

    while abs(sol2[0] - sol2[1]) > tols:

        if helper(log10S2) != (0, 10):

            sol2 = helper(log10S2) 
            log10Shelp = log10S2
            log10S2 = (log10S1 + log10S2) / 2
            print "%s" %str(log10S2)

        else:

            log10S1new = log10S2
            log10S2 = (log10S1 + log10Shelp) / 2
            log10S1 = log10S1new
            
            sol2 = helper(log10S2) 


    return log10Shelp

#--------------------------------------------------------------------------


            
    
