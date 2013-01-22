"""
This script makes the plots makes the "energy_contrib.pdf" and "virial_terms.pdf" plots. First plot shows the contribution of
the various terms in the cooling equation to the global energy balance at the convective boundary & Hill radius. Second plot
shows the terms in the virial equilibrium equation at the convective boundary & Hill radius.
"""

import os, sys
lib_path = os.path.abspath('../')
sys.path.append(lib_path)

import numpy
from numpy import pi
import wr_selfgrav_poly as wr
from userpath import userpath
from disk_model import Td, Pd, rhod
from constants import Me, Re, G, RH, RHe
from parameters import rc, Y, gamma, Cv, delad
import matplotlib.pyplot as plt
from atmprofiles_selfgrav import lmplot, lcontrib, lcontribRH

prof = wr.atmprofileread(userpath + 'dat/MODELS/TEST/test_RH_prof.txt')
param = wr.atmparamread(userpath + 'dat/MODELS/TEST/test_RH_par.txt')

xi = 3 * (gamma - 1)

M = param[:, 0]
Egcb = param[:, 9]
Etotcb = param[:, 11]
Egtot = param[:, 12]
Etot = param[:, 14]
Pcb = param[:, 7]
Tcb = param[:, 8]
Pc = param[:, 5]
rcb = param[:, 2] * Re
RB = param[:, 3] * Re
RHill = param[:, 4] * Re

surfcb = (4 * pi / xi) * (rcb**3 * Pcb - rc**3 * Pc)
surfcbout = (4 * pi / xi) * (rcb**3 * Pcb)
surftot = (4 * pi / xi) * (RHill**3 * Pd - rc**3 * Pc)
surftotout = (4 * pi / xi) * (RHill**3 * Pd)
surfcore = (4 * pi / xi) * (rc**3 * Pc)

    


#plot showing virial terms at radiative-convective boundary and at Hill radius
f, axarr = plt.subplots(2, sharex = True)
axarr[0].semilogx(M, - (xi - 1) / xi * Egcb/10**38, M, -Etotcb/10**38, M, surfcb/10**38, M, surfcbout/10**38)
axarr[0].axhline(y = 0, linestyle = '--', color = 'black')
axarr[0].set_ylabel(r'$E/10^{38}$ (erg)')
axarr[0].legend((r'$-(\xi-1)/\xi E_g$',r'$-E_{\mathrm{tot}}$', r'$(4\pi/\xi)(R_{cb}^3 P_{cb}-R_c^3 P_c)$',\
                 r'$(4\pi/\xi) R_{cb}^3 P_{cb}$'),
                loc = 2, prop = {'size' : 7})
axarr[0].set_title(r'Virial terms at convective boundary: $M_c=5 M_{\oplus}$, $a=10$ AU')
axarr[0].set_xlim(5, 400)
axarr[0].set_ylim(-20, 30)


axarr[1].semilogx(M, - (xi - 1) / xi * Egtot/10**38, M, -Etot/10**38, M, surftot/10**38, M, surftotout/10**38)
axarr[1].set_xlabel(r'$M(M_{\oplus})$')
axarr[1].axhline(y = 0, linestyle = '--', color = 'black')
axarr[1].set_ylabel(r'$E/10^{38}$ (erg)')
axarr[1].legend((r'$-(\xi-1)/\xi E_g$', r'$-E_{\mathrm{tot}}$', r'$(4\pi/\xi)(R_H^3 P_d-R_c^3 P_c)$', \
                r'$(4\pi/\xi) R_{H}^3 P_{d}$'),
                loc = 2, prop = {'size' : 7})
axarr[1].set_title(r'Virial terms at Hill radius: $M_c=5 M_{\oplus}$, $a=10$ AU')
axarr[1].set_ylim(-10, 40)

plt.savefig(userpath + 'figs/ModelAtmospheres/RadSelfGravPoly/RHill_boundary/virial_terms.pdf')
plt.show()


#plot showing the contribution of various energy terms to luminosity

lcb = lcontrib(param[1:], prof[1:])
deltaecb = l[0]
eacccb = l[1]
Pdvcb = l[2]
extraegcb = l[3]

ltot = lcontribRH(param[1:], prof[1:])
deltaetot = ltot[0]
eacctot = ltot[1]
Pdvtot = ltot[2]
extraegtot = ltot[3]

f2, vir = plt.subplots(2, sharex = True)
vir[0].semilogx(M[1:-1], deltaecb / 10**38, M[1:-1], eacccb / 10**38, M[1:-1], Pdvcb / 10**38, M[1:-1], extraegcb / 10**38, \
             M[1:-1], (deltaecb + eacccb + Pdvcb + extraegcb) / 10**38)
vir[1].semilogx(M[1:-1], deltaetot / 10**38, M[1:-1], eacctot / 10**38, M[1:-1], Pdvtot / 10**38, M[1:-1], extraegtot / 10**38, \
             M[1:-1], (deltaetot + eacctot + Pdvtot + extraegtot) / 10**38)
vir[0].axhline(y = 0, linestyle = '--', color = 'black')
vir[1].axhline(y = 0, linestyle = '--', color = 'black')
vir[1].set_xlim(5, 1000)
vir[1].set_ylim(-1, 2)
vir[1].set_xlabel(r'$M (M_{\oplus})$')
vir[0].set_ylabel(r'$E/10^{38}$ (erg)')
vir[1].set_ylabel(r'$E/10^{38}$ (erg)')
vir[0].set_title(r'Cooling terms at the $R_{cb}$: $r=10$ AU and $M_c=5 M_{\oplus}$')
vir[1].set_title(r'Cooling terms at the $R_H$: $r=10$ AU and $M_c=5 M_{\oplus}$')
vir[1].legend((r'$- \Delta{E}$', r'$e_{acc} \Delta{M}$', r'$ - P \Delta{V}$',\
            r'$- \frac{G M}{3 R} \Delta{M}$', r'Total cooling $Ldt$'), prop = {'size' : 8}, loc = 9)
vir[0].legend((r'$- \Delta{E}$', r'$e_{acc} \Delta{M}$', r'$ - P \Delta{V}$',\
            r'$- \frac{G M}{3 R} \Delta{M}$', r'Total cooling $Ldt$'), prop = {'size' : 8}, loc = 9)
plt.savefig(userpath + 'figs/ModelAtmospheres/RadSelfGravPoly/RHill_boundary/energy_contrib.pdf')
plt.show()






