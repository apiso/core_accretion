
"""This module calculates the temperature and pressure corrections at the
convective boundary. The functions are called by shooting_self_grav_mixt_lin."""

"""AY: this docstring needs updating and more detail.  An essential detail is that
it's for use in non-self-grav radiative zones."""


import numpy
import scipy
from scipy import optimize
from interpolation_functions_logt import interpolationp2, interpolationdelad
from table_size import Nt, Np, Ncol
from disk_model import Td, Pd
from parameters import beta

#AY: just cut if this isn't in use.  EOS folder has this code (and better with unpacking instead of [0], [1], etc.)
#EOS=EOSread(Nt, Np, Ncol, '../dat/H_HE_TAB_I_30_EXT.txt')[0]
#T=EOSread(Nt, Np, Ncol, '../dat/H_HE_TAB_I_30_EXT.txt')[1]
#EOSs=EOSread(Nt, 25, 6, '../dat/EOS_S_H_HE_30.txt')[0]

def p(t, S):
    """Interpolates the pressure at temperature t and entropy S"""
    return 10**interpolationp2(numpy.log10(t), S)

def delad(t, S):
    """Interpolates delad at temperature t and entropy S"""
    return interpolationdelad(numpy.log10(t), S)

delinf = 1. / (4 - beta)


def cbcorr(S):
    """
    Calculates temperature and pressure at convective boundary, at entropy S

    output: Tb, Pb (T and P at convective boundary)
    """
    def f(t):
        """Temperature dependence on pressure in the radiative zone"""
        return 1 - (Td/t)**(4 - beta) - (delad(t, S) / delinf) * ( \
            1 - (Pd / p(t, S)))
    Tb = optimize.brentq(f, Td, 5 * Td) #solves for Tb=f(Tb) at the convective boundary
    Pb = p(Tb, S)
    return Tb, Pb



