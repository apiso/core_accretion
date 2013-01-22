"""
Module to test whether structure solutions satisfy virial equilibrium.
See also hbinteract for interactive code to call/use the functions in this module.

So far this is for polytropes that integrate outwards from a core of specified mass and radius.
"""

import numpy as np
from scipy.integrate import odeint
from astropysics.constants import G, kb, mp, Rb, Me, Re
from scipy.integrate import odeint

p4 = 4*np.pi
mu = 2*mp
Rg = kb/mu #gas constant

def DERIVSlnP_liny(y , x , K , gamma = 7./5):
    """
    structure equations dy/yx with
    x = lnP as independent var
    y = (r, M , Eg - Eo, Iu + Eo)
    Eg & Iu -- gravitational energy and Iu (the 3 p/rho dm integral in virial eqn) are diagnostic
    the constant Eo is set by the initial conditions (typically to G*Mc^2 / Rc)
    """

    return np.array([- K**(1 / gamma) / G * y[0]**2 / y[1] * np.exp( x * (1 - 1/gamma)) , \
                     - p4 / G * y[0]**4 / y[1] * np.exp(x) , \
                     p4 * y[0]**3 * np.exp(x) , \
                     - 3 * p4 * K**(1/gamma) / G * y[0]**4 / y[1] * np.exp(x * (2 - 1/gamma))])

def goout(Pc , Pout , Tc , ngrid = 1000. , gamma = 7./5, McE  = 5., rhoco = 3.):
    """
    integrate structure equations along adiabat between two pressure values

    > x , y , K  = goout(... )

    Parameters:
    (Pc , Pout , Tc) -- central, outer pressure and central temp. in cgs

    Keyword parameters:
    ngrid -- points in atmospheric grid
    gamma -- adiabatic index (e.g. 5/3 , 7/5)
    (McE , rhoco) -- core mass in Earth masses , core density in cgs

    Returns:
    x , y , K ...
    """

    #constants

    rhoc = Pc / (Rg * Tc) #central density
    K = Pc / rhoc**gamma # constant for polytrope

    M0 = McE * Me
    r0 = (3 * M0 / p4 / rhoco)**(1./3)
    E0 = G * M0**2 / r0
    #y0 = (np.log(r0) , np.log(M0) , -E0 , E0)
    y0 = (r0 , M0 , -E0 , E0)

    #set grid and integrate
    x = np.linspace(np.log(Pc) , np.log(Pout) , ngrid)
    y = odeint(DERIVSlnP_liny, y0, x, args = (K, gamma))

    #vircheck(x , y , K)

    return x, y, K

def gooutrho(rhoc , rhoout , fv , ngrid = 1000. , gamma = 7./5, McE  = 5., rhoco = 3.):
    """
    wrapper for goout

    Parameters:
    core and outer density in cgs and core temperature as a fraction of virial value

    Keywords: match go out, could pass all extra arguments with ** trick)
    """
    #constants
    delad = 1 - 1/gamma

    #inital (core) values)
    M0 = McE * Me
    r0 = (3 * M0 / p4 / rhoco)**(1./3)
    Tc = fv * delad * G * M0 / r0 / Rg #core temp as fract of virial
    Pc = rhoc * Rg * Tc

    K = Pc / rhoc**gamma # constant for polytrope
    Pout = K * rhoout**gamma #outer pressure

    return goout(Pc , Pout ,  Tc , ngrid , gamma , McE , rhoco)

def vircheck(x , y):
    """
    check viral balance of an atmospheric solution
    
    Parameters:
    x = ln(P) (in cgs)
    y = (r, m, Eg - E0, Iu - E0) -- Iu being 3P/rho integral in virial thm.
    """

    #E0 = G * np.exp(2 * y[0,1] - y[0,0])
    #Eg = -(np.exp(y[-1,2]) - E0)
    #Iu = np.exp(y[-1,3]) - E0
    E0 = G * y[0,1]**2 / y[0,0]
    Eg = y[-1,2] + E0
    Iu = y[-1,3] - E0
    
    #bcore = p4 * np.exp(3 * y[0 , 0] + x[0])
    #bsurf = -p4 * np.exp(3 * y[-1 , 0] + x[-1])
    bcore = p4 * y[0 , 0]**3 * np.exp(x[0])
    bsurf = -p4 * y[-1 , 0]**3 * np.exp(x[-1])
    
    bal = (Eg + Iu + bcore + bsurf) / Eg
    terms = (-Iu/Eg , -bcore/Eg , -bsurf/Eg)

    print "virial balance within" , bal
    print "ratio of terms to -Eg" , terms

    return

#these integrators were improved on...

def DERIVSlnP(y , x , K , gamma = 7./5):
    """
    structure equations dy/yx with
    x = lnP as independent var
    y = (lnr, lnM)
    Eg & Iu -- gravitational energy and Iu (the 3 p/rho dm integral in virial eqn) are diagnostic
    """

    return np.array([- K**(1 / gamma) / G * np.exp( x * (1 - 1/gamma) + y[0] - y[1]) , \
                     - p4 / G * np.exp(x + 4*y[0] - 2 * y[1]) ])

def DERIVSlnP_vc(y , x , K , gamma = 7./5):
    """
    structure equations dy/yx with
    x = lnP as independent var
    y = (lnr, lnM , Eg - Eo, Iu + Eo)
    Eg & Iu -- gravitational energy and Iu (the 3 p/rho dm integral in virial eqn) are diagnostic
    the constant Eo is set by the initial conditions (typically to G*Mc^2 / Rc)
    """

    return np.array([- K**(1 / gamma) / G * np.exp( x * (1 - 1/gamma) + y[0] - y[1]) , \
                     - p4 / G * np.exp(x + 4*y[0] - 2 * y[1]) , \
                     p4 * np.exp(x + 3*y[0]) , \
                     - 3 * p4 * K**(1/gamma) / G * np.exp(x * (2 - 1/gamma) + 4 * y[0] - y[1])])
