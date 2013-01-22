import os, sys
lib_path = os.path.abspath('../')
sys.path.append(lib_path)

import numpy
import wr_selfgrav_poly as wr
from userpath import userpath
from disk_model import Td, Pd, rhod
from constants import Me, Re
from parameters import rc
import matplotlib.pyplot as plt
from atmprofiles_selfgrav import lmplot

prof = wr.atmprofileread(userpath + 'dat/MODELS/TEST/test_RH_prof.txt')
param = wr.atmparamread(userpath + 'dat/MODELS/TEST/test_RH_par.txt')

tpro = prof[4]
tpar = param[4]

r = tpro[:,0] / Re
P = tpro[:,1]
T = tpro[:,2]
rho = tpro[:,3]
m = tpro[:,4]

M = param[:,0]
#M = tpar[0]
rcb = tpar[2]
RB = tpar[3]
RHill = tpar[4]

#plt.figure(1)
f, axarr = plt.subplots(3, sharex=True)
ax = axarr[2].xaxis

axarr[0].loglog(r, (P - Pd) / Pd)
axarr[0].set_ylabel(r'($P$ - $P_d$)/$P_d$')
axarr[0].set_title(r'$M_c=5.0 M_{\oplus}$, $M_{\mathrm{tot}}=8.06 M_{\oplus}$, $a=10.0$ AU')
axarr[0].yaxis.set_ticks([10**(-2), 10, 10**4, 10**7, 10**10])

axarr[0].axvline(x = rc / Re, linestyle = ':', color = 'black')
axarr[0].axvline(x = rcb, linestyle = '-.', color = 'red')
axarr[0].axvline(x = RB, linestyle = '--', color = 'green')
axarr[0].axvline(x = RHill, linestyle = ':', color = 'black')

axarr[1].loglog(r, (rho - rhod) / rhod)
axarr[1].set_ylabel(r'($\rho$ - $\rho_d$)/$\rho_d$')
axarr[1].yaxis.set_ticks([10**(-2), 10, 10**4, 10**7])

axarr[1].axvline(x = rc / Re, linestyle = ':', color = 'black')
axarr[1].axvline(x = rcb, linestyle = '-.', color = 'red')
axarr[1].axvline(x = RB, linestyle = '--', color = 'green')
axarr[1].axvline(x = RHill, linestyle = ':', color = 'black')

axarr[2].loglog(r[:-1], (T[:-1] - Td) / Td)
axarr[2].set_ylabel(r'($T$ - $T_d$)/$T_d$')
axarr[2].set_xlabel(r'$r(R_{\oplus})$')
axarr[2].yaxis.set_ticks([10**(-6), 10**(-3), 1, 10**3])

axarr[2].axvline(x = rc / Re, linestyle = ':', color = 'black')
axarr[2].axvline(x = rcb, linestyle = '-.', color = 'red')
axarr[2].axvline(x = RB, linestyle = '--', color = 'green')
axarr[2].axvline(x = RHill, linestyle = ':', color = 'black')

axarr[2].annotate(r'$R_c$', xy = (rc / Re, 10**(-3)))
axarr[2].annotate(r'$R_{\mathrm{cb}}$', xy = (rcb, 10**(-3)))
axarr[2].annotate(r'$R_B$', xy = (RB, 10**(-3)))
axarr[2].annotate(r'$R_H$', xy = (RHill, 10**(-3)))

plt.savefig(userpath + 'figs/model atmospheres/RadSelfGravPoly/RHill_boundary/PTrho_vs_r.pdf')
plt.show()


lmplot(param[:145])
plt.axvline(x = M[18], linestyle = '--', color = 'black')
plt.axvline(x = M[36], linestyle = '--', color = 'green')
plt.axvline(x = M[123], linestyle = '--', color = 'red')
plt.annotate(r'$R_B<R_H<H$', rotation = 90, color = 'black', xy = (10, 5.3 * 10**26))
plt.annotate(r'$R_H<R_B<H$', rotation = 90, color = 'black', xy = (25, 5.3 * 10**26))
plt.annotate(r'$R_H<H<R_B$', rotation = 90, color = 'black', xy = (63, 5.3 * 10**26))
plt.annotate(r'$H<R_H<R_B$', rotation = 90, color = 'black', xy = (105, 5.3 * 10**26))
plt.savefig(userpath + 'figs/model atmospheres/RadSelfGravPoly/RHill_boundary/LMplot_v2.pdf')
plt.show()

