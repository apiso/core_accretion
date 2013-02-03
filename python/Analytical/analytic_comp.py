"""This module computes analytic estimates for the core accretion model; the
abalytic model does not take into account self-gravity in either the convective
or radiative region.
Useful for comparison with the full numerical model."""

import numpy as np
from utils.constants import G, kb, mp, Rb, Me, Re, Msun, RH, RHe, sigma, \
     cmperau, RHill, gammafn, mufn, Rfn, Cvfn, kdust, Tdisk, Pdisk, params
from utils.parameters import mstar, FSigma, FT
from scipy.integrate import odeint

def analytic_sols(delad, Y, rhoc, a, Mc, Pcb, beta = 2, theta = 0.556):

    """
    Computes analytic estimates for a given set of gas, disk and core conditions
    and a given value (or set of values) at the radiative convective boundary.

    WARNING: Though this function takes delad as a parameter, the results for
    cooling time and M_CB assume delad = 2./7.  This should be changed.

    Input
    -----
    delad, Y:
        adiabatic index and He mass fraction
    rhoc:
        core density in g/cm^3
    a:
        semi-major axis in AU
    Mc:
        core mass in g
    Pcb:
        pressure at the RCB in dyne/cm^2; can be a scalar or an array
    beta:
        temperature power law index for sut opacity; default 2
    theta:
        correction factor for the pressure at the RCB; default 0.556 for beta=2

    Output
    ------
    Tcb:
        temperature at the RCB
    Matm:
        atmosphere mass (non self-gravitating) in g
    Lcb:
        luminosity at the RCB in erg/s
    tcool:
        cooling time in years
        
    
    """

    p4 = 4 * np.pi
    delinf = 1. / (4 - beta)

    rc = (3 * Mc / (4 * np.pi * rhoc))**(1./3)
    
    R = Rfn(Y)

    Td = Tdisk(a, FT)
    Pd = Pdisk(a, mstar, FSigma, FT)

    chi = (1 - delad / delinf)**(-1. / (4 - beta))
    Tcb = Td * chi
    RB = G * Mc / (R * Td)
    RBprime = (delad / chi) * RB
    PM = (4 * delad**(3./2) / (5 * np.pi**2 * np.sqrt(chi))) * (\
        G * Mc**2 / RBprime**4)

    Mcb = Mc * (Pcb / PM) / (np.sqrt(np.log(Pcb / (theta * Pd))))

    """While the true L at convective zone depends on Mcb, the analytic
    theory only considers Mc to make integrals easier"""
    ML = Mc #+ Mcb neglected
    Lo = 64 * np.pi * G * ML * sigma * Td**4 /\
               (3 * kdust(Td) * Pd) * delad * chi**(4 - beta)

    Lcb = Lo * Pd / Pcb

    tcool = 4 * np.pi * Pcb**2 / Pd * RBprime**(7./2) / (Lo * np.sqrt(rc))

    #Now add correction to atmospheric mass for radiative zone mass.
    #Alternately the comparison of numerical and analytic could be done at
    #RCB

    rcb = RB / np.log(Pcb/(theta*Pd))
    """Below is only a weakly self-gravitating Hill radius as the radiative
    zone mass is excluded, and further `Mcb` is from a non-self-grav. calc.
    Good enough for our purposes here."""
    rHill = RHill(Mc + Mcb, a, Ms=mstar*Msun)
         
    def mPstruct(y, r, Tfunc):
         """HB struct equations requiring `Tfunc` = T(P).
         Putting as sub-fn for now to avoid passing params."""

         m, P = y
         T = Tfunc(P)
         rho = P / R / T
         dmdr = p4 * rho * r**2
         dPdr = - rho * G * m / r**2
         return dmdr, dPdr

    def TPrad(P):
         """T-P relation that requires setting Pcb value outside fn
         and knowing Pd.  Classes are DEFINITELY the better way to do
         this, but we're not there yet."""

         myPcb = TPrad.Pcb
         #WARNING: in this form the temp is specific to alpha = 0,
         #         no pressure dependence in powerlaw opacity
         Tfac = (1 + (P - Pd) / myPcb / (delinf/delad - 1))**delinf
         return Td * Tfac

    #if only a single `Pcb` is passed make quantities iterable
    if not hasattr(Pcb, '__iter__'):
          Pcb = [Pcb]
          rcb, Mcb, rHill = [rcb], [Mcb], [rHill]

    Mrad = np.empty_like(Pcb) 
    for i, Pcb_i in enumerate(Pcb):
         TPrad.Pcb = Pcb_i
         y = odeint(mPstruct, (Mc + Mcb[i], Pcb_i), (rcb[i], RB, rHill[i]),
                    args=(TPrad,))
         Mtot = y[1, 0] #for row entry (first) use 1 for RB mass
                        #and 2 (or -1) for Hill
         Mrad[i] = Mtot - (Mc + Mcb[i])
    
    return Tcb, Mcb, Lcb, tcool / (365 * 24 * 3600), rcb, Mrad, rHill
    #if desired could switch Mcb to total Matm above.






