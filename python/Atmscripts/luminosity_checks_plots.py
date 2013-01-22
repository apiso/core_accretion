"""
This script makes the 'neglected_luminsoty.pdf' plot showing how the estimated dL/dm varies
throughout the radiative region.
"""

import os, sys
#lib_path = os.path.abspath('../')
#sys.path.append(lib_path)

import numpy
from numpy import pi
import wr_selfgrav_poly as wr
#from userpath import userpath
from disk_model import Td, Pd, rhod
from constants import Me, Re, G, RH, RHe
from parameters import rc, Y, gamma, Cv, delad
import matplotlib.pyplot as plt
from atmprofiles_selfgrav import luminositycheck

prof = wr.atmprofileread(userpath + 'dat/MODELS/TEST/test_RH_prof.txt')
param = wr.atmparamread(userpath + 'dat/MODELS/TEST/test_RH_par.txt')

m = luminositycheck(prof[10:995], param[10:995])[-1]
dLdm = luminositycheck(prof[10:995], param[10:995])[-2]

vec = [10, 15, 20]
Mtot = [param[20,0], param[25, 0], param[30, 0]]

plt.figure()
for i in vec:
    plt.plot(m[i] / Me, dLdm[i])
plt.legend((r'$M_{\mathrm{tot}}= % s M_{\oplus}$' % str(Mtot[0]), \
            r'$M_{\mathrm{tot}}= % s M_{\oplus}$' % str(Mtot[1]), \
             r'$M_{\mathrm{tot}}= % s M_{\oplus}$' % str(Mtot[2])))
plt.axhline(y = 0, linestyle = '--', color = 'black')
plt.xlabel(r'$m (M_{\oplus})$')
plt.ylabel(r'$dL/dm$ (erg s$^{-1}$ g$^{-1}$)')
plt.title(r'$dL/dm$ between $CB$ and $R_H$: $a=10$ AU, $M_c=5M_{\oplus}$')
#plt.savefig(userpath + 'figs/ModelAtmospheres/RadSelfGravPoly/RHill_boundary/neglected_luminosity.pdf')
plt.show()




