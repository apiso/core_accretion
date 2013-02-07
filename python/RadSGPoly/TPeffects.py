"""

Generates series of atmosphere profiles for a set of gas, disk and core
conditions and writes them in .npz files. Also contains functions to load
profiles.

"""



from utils.constants import G, kb, mp, Rb, Me, Re, Msun, RH, RHe, sigma, \
     cmperau, RHill, gammafn, mufn, Rfn, Cvfn, kdust, Tdisk, Pdisk, params
from utils.parameters import mstar, FSigma, FT
from utils.userpath import userpath
import numpy
import scipy
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from types import FunctionType as function
from collections import namedtuple
from scipy import integrate
from utils import cb_corrections as cb
from utils.interpolation_functions_logp import interpolations
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.optimize import fminbound, brentq
#from shooting_poly import Ltop, shoot, prms
from profiles_poly import profiles_write as pp

def Tdfixed(a, delad, Y, rhoc, Mc, fact):

    rc = (3 * Mc / (4 * numpy.pi * rhoc))**(1./3)
    prms = params(Mc, rc, a, delad, Y, gamma = gammafn(delad), R = Rfn(Y), \
                  Cv = Cvfn(Y, delad), Pd = Pdisk(a, mstar, FSigma, FT) * fact, \
                  Td = Tdisk(a, FT), kappa = kdust)

    pp(500, 300, 10**21, 10**30, 1.2*Mc, 2.5*Mc, 'Mc' + \
                     str(Mc/Me)[:4], prms = prms, disk = 0)

def Pdfixed(a, delad, Y, rhoc, Mc, fact):

    rc = (3 * Mc / (4 * numpy.pi * rhoc))**(1./3)
    prms = params(Mc, rc, a, delad, Y, gamma = gammafn(delad), R = Rfn(Y), \
                  Cv = Cvfn(Y, delad), Pd = Pdisk(a, mstar, FSigma, FT), \
                  Td = Tdisk(a, FT) * fact, kappa = kdust)

    pp(500, 300, 10**21, 10**30, 1.2*Mc, 2.5*Mc, 'Mc' + \
                     str(Mc/Me)[:4], prms = prms, disk = 0)
