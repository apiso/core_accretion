os.chdir(userpath + '/python')

import numpy
from utils.constants import Me, Re
import matplotlib.pyplot as plt
from utils.userpath import userpath
from RadSGRealGas.shooting import prms
from RadSGRealGas.profiles import atmload
from RadSGRealGas.cooling import cooling_global as cg, cooling_local as cl

os.chdir(userpath + '/python/Atmscripts')


load = 1
savefig = 0
savefigt = 0
out = 'rcb'

if load != 0:

    atm = atmload('Mc6.0.npz', prms = prms)
    model = atm[0]
    param = atm[1]
    prof = atm[2]

    Mass = param.MB
    #M = M[1:]

    cool = cg(param, prof, prms = prms, out = out)
    t = cool[0]
    deltae = cool[2]
    eaccdm = cool[3]
    PdV = cool[4]
    Lglobdt = -deltae + eaccdm - PdV

    Llocdt = cl(param, prof, prms = prms, out = out)
    Llocdtnegl = cl(param, prof, prms = prms, out = 'RHill', onlyrad = 1)

L = param.L
dM = (Mass[1:] - Mass[:-1]) * Me
M = Mass[:-1] - model.Mco / Me #only atmosphere mass for plot

plt.loglog(M, - deltae / dM, 'r', label = r'$-dE/dM$')
plt.loglog(M, deltae / dM, '--r')
plt.loglog(M, eaccdm / dM, 'g', label = r'$-e_{\mathrm{acc}}$')
plt.loglog(M, - eaccdm / dM, '--g')
plt.loglog(M, - PdV / dM, 'y', label = r'$-PdV/dM$')
plt.loglog(M, PdV / dM, '--y')
plt.loglog(M, Lglobdt / dM, 'black', label = r'$L_{\mathrm{glob}}dt/dM$')
plt.loglog(M, - Lglobdt / dM, 'black', linestyle = '--')
plt.loglog(M, Llocdt / dM, 'b', label = r'$L_{\mathrm{loc}}dt/dM$')
plt.loglog(M, - Llocdt / dM, 'b', linestyle = '--')
plt.loglog(M, Llocdtnegl / dM, 'cyan', label = r'$L_{\mathrm{negl}}dt/dM$')
plt.loglog(M, - Llocdtnegl / dM, 'cyan', linestyle = '--')
plt.loglog(M, L[:-1] / 10**17, 'orange', label = r'$L/10^{17}$')
plt.legend(loc = 3, prop = {'size':10})
plt.xlabel(r'$M_{\mathrm{atm}}(M_{\oplus})$')
#plt.xlim(0.1,1000)
if savefig != 0:
    plt.savefig(userpath + '/figs/ModelAtmospheres/RadSelfGravPoly/cooling_a10_Mc5_' + out + '.pdf')
plt.show()



tcb = cg(param, prof, prms = prms, out = 'rcb')[0] / (365 * 24 * 3600)
tHill = cg(param, prof, prms = prms, out = 'RHill')[0] / (365 * 24 * 3600)
tB = cg(param, prof, prms = prms, out = 'RB')[0] / (365 * 24 * 3600)


n = len(tcb)
tcumcb = 0 * numpy.ndarray(shape = (n), dtype = float)
tcumHill = 0 * numpy.ndarray(shape = (n), dtype = float)
tcumB = 0 * numpy.ndarray(shape = (n), dtype = float)

for i in range(n):
    tcumcb[i] = sum(tcb[:i])
    tcumHill[i] = sum(tHill[:i])
    tcumB[i] = sum(tB[:i])
    
plt.loglog(M, tcb, 'b', label = 'CB')
plt.loglog(M, tHill, 'g', label = 'Hill')
plt.loglog(M, tB, 'r', label = 'Bondi')
plt.loglog(M, -tcb, '--b')
plt.loglog(M,- tHill, '--g')
plt.loglog(M,- tB, '--r')
plt.xlabel(r'$M_{\mathrm{atm}}(M_{\oplus})$')
plt.ylabel(r'$\Delta t$ (yrs)')
plt.legend()
if savefigt != 0:
    plt.savefig(userpath + '/figs/ModelAtmospheres/RadSelfGravPoly/cooling_time_a10_Mc5.pdf')
plt.show()


plt.loglog(tcumcb, M, 'b', label = 'CB')
plt.loglog(tcumHill, M, 'g', label = 'Hill')
plt.loglog(tcumB, M, 'r', label = 'Bondi')
plt.loglog(-tcumcb, M, '--b')
plt.loglog(- tcumHill, M, '--g')
plt.loglog(- tcumB, M, '--r')
plt.ylabel(r'$M_{\mathrm{atm}}(M_{\oplus})$')
plt.xlabel(r'$t$ (yrs)')
plt.title('Cumulative cooling time')
plt.legend()
if savefigt != 0:
    plt.savefig(userpath + '/figs/ModelAtmospheres/RadSelfGravPoly/cumulative_cooling_time_a10_Mc5.pdf')
plt.show()


