"""This module computes analytic estimates for the core accretion model; the
abalytic model does not take into account self-gravity in either the convective
or radiative region.
Useful for comparison with the full numerical model."""

import numpy as np
from utils.constants import G, kb, mp, Rb, Me, Re, Msun, RH, RHe, sigma, \
     cmperau, RHill, gammafn, mufn, Rfn, Cvfn, kdust, Tdisk, Pdisk, params
from utils.parameters import mstar, FSigma, FT

def analytic_sols(delad, Y, rhoc, a, Mc, Pcb, beta = 2, theta = 0.556):

    """
    Computes analytic estimates for a given set of gas, disk and core conditions
    and a given value (or set of values) at the radiative convective boundary.

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
    
    delinf = 1. / (4 - beta)

    rc = (3 * Mc / (4 * np.pi * rhoc))**(1./3)

    R = Rfn(Y)

    Td = Tdisk(a, FT)
    Pd = Pdisk(a, mstar, FSigma, FT)

    chi = (1 - delad / delinf)**(-1. / (4 - beta))
    Tcb = Td * chi
    RBprime = (delad / chi) * (G * Mc / (R * Td))
    PM = (4 * delad**(3./2) / (5 * np.pi**2 * np.sqrt(chi))) * (\
        G * Mc**2 / RBprime**4)

    Matm = Mc * (Pcb / PM) / (np.sqrt(np.log(Pcb / (theta * Pd))))

    Lo = 64 * np.pi * G * (Mc + Matm) * sigma * Td**4 /\
               (3 * kdust(Td) * Pd) * delad * chi**(4 - beta)

    Lcb = Lo * Pd / Pcb

    tcool = 4 * np.pi * Pcb**2 / Pd * RBprime**(7./2) / (Lo * np.sqrt(rc))

    return Tcb, Matm, Lcb, tcool / (365 * 24 * 3600)







