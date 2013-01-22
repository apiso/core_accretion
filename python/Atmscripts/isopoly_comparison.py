import numpy
from IsoPoly import isopoly as ip
from IsoPoly.energetics import virialterms
from utils.userpath import userpath
from utils.constants import Me, Re
from utils.parameters import Mc
from RadNoSGPoly.profiles_poly import paramload, profload
#from RadNoSGPoly.atmprofiles_poly import prof_merge, param_merge
import matplotlib.pyplot as plt

load = 0

#arr = ip.sequence(nP = 500, mco = 5*Me)
#vir = virialterms(file = "mc5au10_N1000.npz")

##paramlow = wr.atmparamread(userpath + '/dat/MODELS/RadNoSGPoly/Mc5_low_param_iso.txt')[0]
##proflow = wr.atmprofileread(userpath + '/dat/MODELS/RadNoSGPoly/Mc5_low_prof_iso.txt')[0]
##paramhigh = wr.atmparamread(userpath + '/dat/MODELS/RadNoSGPoly/Mc5_high_param_iso.txt')[0]
##profhigh = wr.atmprofileread(userpath + '/dat/MODELS/RadNoSGPoly/Mc5_high_prof_iso.txt')[0]
##
##param = param_merge(paramlow, paramhigh)
##prof = prof_merge(proflow, profhigh)

if load != 0:

    param = paramload('Mc5param_logp_iso_500.npz')
    prof = profload('Mc5prof_logp_iso_500.npz')
    npzdat = numpy.load(userpath + '/dat/MODELS/isopoly/Hillmc5au10.npz')
    arr = npzdat['atmser'].view(numpy.recarray)
    vir = virialterms('Hillmc5au10.npz')

Mcb = param.Mcb
Pcb = param.Pcb
Pc = param.Pc
Rcb = param.rcb
Eg = param.Eg
Etot = param.Etot

dEvir = vir.Etot[1:] - vir.Etot[:-1]
dE = Etot[1:] - Etot[:-1]

dMcb = (Mcb[1:]-Mcb[:-1]) * Me
dMcbarr = (arr.Mcb[1:] - arr.Mcb[:-1]) 

Mcb = Mcb - Mc / Me
Mcbarr = (arr.Mcb - Mc) / Me

f, axarr = plt.subplots(2, sharex = True)
axarr[0].loglog(Mcb, Pcb, Mcbarr, arr.Pcb, 'r')
axarr[0].set_ylabel(r'$P_{CB}$ (dyn cm$^{-2}$)')
axarr[0].legend(('AP', 'AY'))
axarr[1].loglog(Mcb, Rcb, Mcbarr, arr.Rcb / Re, 'r')
axarr[1].set_ylabel(r'$R_{CB}$ ($R_{\oplus}$)')
axarr[1].set_xlabel(r'$M_{\mathrm{atm}} (M_{\oplus})$')
plt.savefig(userpath + '/figs/ModelAtmospheres/isopoly/compare_struct_Hill.pdf')
plt.show()

f2, ax2 = plt.subplots(2, sharex = True)
ax2[0].loglog(Mcb, - Etot)
ax2[0].loglog(Mcbarr, - vir.Etot, 'r')
ax2[0].loglog(Mcb, Etot, '--b')
ax2[0].loglog(Mcbarr, vir.Etot, '--r')
ax2[0].set_ylabel(r'$-E_{tot}$ (erg)')
ax2[0].legend(('AP', 'AY'), loc = 2)

ax2[1].loglog(Mcb[:-1], - dE / dMcb, Mcb[:-1], dE / dMcb, '--b')
ax2[1].loglog(Mcbarr[:-1], - dEvir / dMcbarr, 'r')
ax2[1].loglog(Mcbarr[:-1], dEvir / dMcbarr, '--r')
ax2[1].set_ylabel(r'$-dE/dM$ (erg g$^{-1}$)')
ax2[1].set_xlabel(r'$M_{\mathrm{atm}} (M_{\oplus})$')
plt.savefig(userpath + '/figs/ModelAtmospheres/isopoly/compare_energy_Hill.pdf')

plt.show()







