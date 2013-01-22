
"""
This module calculates the temperature and pressure in the disk, based on the
disk model of Chiang & Youdin (2010). It also calculated temperature and pressure
corrections at the convective boundary, with the assumption that Pb >> Pd
"""

import numpy
import scipy
from scipy import optimize
from EOS import EOSread
from constants import Tdisk, Pdisk
from parameters import a, FT, FSigma, mstar, beta, R
from table_size import Nt, Np, Ncol
from userpath import userpath
from interpolation_functions_logp import interpolations

EOS=EOSread.EOSread(Nt, Np, Ncol, userpath + '/dat/EOStables/H_HE_TAB_I_30_EXT.txt')[0]
T=EOSread.EOSread(Nt, Np, Ncol, userpath + '/dat/EOStables/H_HE_TAB_I_30_EXT.txt')[1]
EOSs=EOSread.EOSread(Nt, 25, 6, userpath + '/dat/EOStables/EOS_S_H_HE_30.txt')[0]

def k(T, b):
    """Opacity (generally applicable for beta = 2)"""
    return 2 * (T/100)**b

def kappa(T):
    """Opacity for a given beta from the parameters file (default beta = 2)"""
    return k(T, beta)

#Td = 30.0
Td = Tdisk(a, FT) #Td for given r and Ft from parameters; imported in shooting code
#Pd = 0.6
Pd = Pdisk(a, mstar, FSigma, FT) #Pd for given r, mstar, Fsigma, Ft from parameters
rhod = Pd / (R * Td) #disk density
Sd = interpolations(numpy.log10(Td), numpy.log10(Pd))












