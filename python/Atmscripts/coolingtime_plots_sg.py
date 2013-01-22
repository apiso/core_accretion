import os
from utils.userpath import userpath
os.chdir(userpath + '/python')

import numpy
from utils.constants import Me, Re, gammafn, Rfn, Cvfn, Pdisk, Tdisk, kdust, \
     params
from utils.parameters import FSigma, FT, mstar
import matplotlib.pyplot as plt

from RadSGPoly.profiles_poly import atmload
from RadSGPoly.cooling_poly import cooling_global as cg, cooling_local as cl, critical

os.chdir(userpath + '/python/Atmscripts')

delad = 2./7
Y = 0.3

load = 1
cool = 1
savefig = 0

#a = 10.0

prms = params(Me, Re, a, delad, Y, gamma = gammafn(delad), R = Rfn(Y), \
              Cv = Cvfn(Y, delad), Pd = Pdisk(a, mstar, FSigma, FT), \
              Td = Tdisk(a, FT), kappa = kdust)

if load != 0:

    atm5 = atmload('Mc1.0.npz', prms = prms)
    atm6 = atmload('Mc2.0.npz', prms = prms)
    atm7 = atmload('Mc3.0.npz', prms = prms)
    atm8 = atmload('Mc4.0.npz', prms = prms)
    atm9 = atmload('Mc5.0.npz', prms = prms)
    atm10 = atmload('Mc6.0.npz', prms = prms)
    atm11 = atmload('Mc7.0.npz', prms = prms)
    atm12 = atmload('Mc8.0.npz', prms = prms)
    atm13 = atmload('Mc9.0.npz', prms = prms)
    atm14 = atmload('Mc10.0.npz', prms = prms)
    #atm15 = atmload('Mc13.0.npz', prms = prms)

Mc = numpy.linspace(1, 10, 10)

model = [atm5[0], atm6[0], atm7[0], atm8[0], atm9[0], atm10[0], atm11[0], atm12[0], atm13[0], atm14[0]]
param = [atm5[1], atm6[1], atm7[1], atm8[1], atm9[1], atm10[1], atm11[1], atm12[1], atm13[1], atm14[1]]
prof = [atm5[2], atm6[2], atm7[2], atm8[2], atm9[2], atm10[2], atm11[2], atm12[2], atm13[2], atm14[2]]

n = len(model)

if cool != 0:

    t = 0 * numpy.ndarray(shape = (n), dtype = float)
    
    for i in range(n):
        temp = critical(param[i], prof[i], model[i])
        dt = temp[-1]
        t[i] = sum(dt)
    
plt.plot(Mc, t / (365 * 24 * 3600) / 10**6, '-o')
plt.axhline(y = 3, linestyle = '--', color = 'black', label = 't = 3 Myrs')
plt.xlabel(r'$M_c (M_{\bigoplus})$ ')
plt.ylabel(r'$t_{\mathrm{cool}}$(Myrs)')
plt.legend()
plt.title('a=' + str(int(a)) + ' AU')
if savefig != 0:
    plt.savefig(userpath + '/figs/ModelAtmospheres/RadSelfGravPoly/coolingtime_vs_Mc_' + str(int(a)) + 'au.pdf')
plt.show()



