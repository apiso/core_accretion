import numpy
from utils.constants import Me, Re
import matplotlib.pyplot as plt
from utils.userpath import userpath
from RadSGPoly.shooting_poly import params, prms
from RadSGPoly.profiles_poly import atmload
from RadSGPoly.cooling_poly import critical, cooling_global as cg, \
     cooling_local as cl

load = 0
savefig = 0
savefigt = 0
out = 'rcb'

if load != 0:

    atm = atmload('Mc5.0.npz', prms = prms)
    model = atm[0]
    param = atm[1]
    prof = atm[2]


    crit = critical(param, prof, model)

    param = crit[0]
    prof = crit[1]

    Mass = param.MB
    
    cool = cg(param, prof, model, out = out)

    dt = cool[0]
    deltae = cool[2]
    eaccdm = cool[3]
    PdV = cool[4]
    Lglobdt = - deltae + eaccdm - PdV

    Llocdt = cl(param, prof, prms = prms, out = out)
    Llocdtnegl = cl(param, prof, prms = prms, out = 'RHill', onlyrad = 1)

L = param.L
dM = (Mass[1:] - Mass[:-1]) * Me
M = Mass[:-1] - model.Mco / Me #only atmosphere mass for plot

plt.semilogy(M, - deltae / dt, 'r', label = r'$-dE/dt$')
plt.loglog(M, deltae / dt, '--r')
plt.loglog(M, eaccdm / dt, 'g', label = r'$-e_{\mathrm{acc}} dM/dt$')
plt.loglog(M, - eaccdm / dt, '--g')
plt.loglog(M, - PdV / dt, 'y', label = r'$-PdV/dt$')
plt.loglog(M, PdV / dt, '--y')
#plt.loglog(M, Lglobdt / dt, 'black', label = r'$L_{\mathrm{glob}}$')
#plt.loglog(M, - Lglobdt / dt, 'black', linestyle = '--')
#plt.loglog(M, Llocdt / dt, 'b', label = r'$L_{\mathrm{loc}}$')
#plt.loglog(M, - Llocdt / dt, 'b', linestyle = '--')
plt.loglog(M, Llocdtnegl / dt, 'black', label = r'$L_{\mathrm{negl}}$')
plt.loglog(M, - Llocdtnegl / dt, 'black', linestyle = '--')
plt.loglog(M, L[:-1], 'b', label = r'$L$')
plt.legend(loc = 3)
plt.xlabel(r'$M_{\mathrm{atm}}(M_{\oplus})$')
#plt.xlim(0.1,1000)
if savefig != 0:
    plt.savefig(userpath + '/figs/ModelAtmospheres/RadSelfGravPoly/cooling_a10_Mc5_' + out + '.pdf')
plt.show()




