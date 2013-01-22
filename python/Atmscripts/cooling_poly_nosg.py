import numpy
from utils.parameters import Mc
from utils.constants import Me, Re
import matplotlib.pyplot as plt
from utils.userpath import userpath
from RadNoSGPoly.profiles_poly import paramload, profload
from RadNoSGPoly.cooling_poly import cooling_global as cg, cooling_local as cl

load = 1
savefig = 0

if load != 0:

    param = paramload('Mc5param_logp_2000.npz')
    prof = profload('Mc5prof_logp_2000.npz')

    M = param.Mcb
    #M = M[1:]

    cool = cg(param, prof)
    t = cool[0]
    deltae = cool[2]
    eaccdm = cool[3]
    PdV = cool[4]
    Lglobdt = -deltae + eaccdm - PdV

    Llocdt = cl(param, prof)

dM = (M[1:] - M[:-1]) * Me
M = M[:-1] - Mc / Me #only atmosphere mass for plot

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
plt.legend(loc = 1)
#plt.xlim(0.2,0.3)
if savefig != 0:
    plt.savefig(userpath + '/figs/ModelAtmospheres/Rad_Compare/cooling_a10_Mc5_5000.pdf')
plt.show()


