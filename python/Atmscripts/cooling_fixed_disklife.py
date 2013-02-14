import os
os.chdir('/home/apiso/repos/core_accretion/python')
#os.chdir('/Users/ana-mariapiso/core_accretion/python')

import numpy
from utils.constants import Me, Re, gammafn, Rfn, Cvfn, Pdisk, Tdisk, kdust, \
     params, mufn
from utils.userpath import userpath
from utils.parameters import FSigma, FT, mstar
import matplotlib.pyplot as plt
from RadSGPoly import plot_helper as ph
from RadSGPoly.profiles_poly import atmload
import matplotlib.cm as cm

a = 5.0
delad27 = 2./7
delad25 = 2./5
Y00 = 0.0
Y03 = 0.3

prms00 = params(Me, Re, a, delad27, Y00, gamma = gammafn(delad27), R = Rfn(Y00), \
              Cv = Cvfn(Y00, delad27), Pd = Pdisk(a, mstar, FSigma, FT), \
              Td = Tdisk(a, FT), kappa = kdust)

prms03 = params(Me, Re, a, delad27, Y03, gamma = gammafn(delad27), R = Rfn(Y03), \
              Cv = Cvfn(Y03, delad27), Pd = Pdisk(a, mstar, FSigma, FT), \
              Td = Tdisk(a, FT), kappa = kdust)

prms03m = params(Me, Re, a, delad25, Y03, gamma = gammafn(delad25), R = Rfn(Y03), \
              Cv = Cvfn(Y03, delad25), Pd = Pdisk(a, mstar, FSigma, FT), \
              Td = Tdisk(a, FT), kappa = kdust)




Mct_a10_00 = ph.t_vs_Mc_fixed_a(2. / 7, 0.0, a, returnt = 1)
Mct_a10_03 = ph.t_vs_Mc_fixed_a(2. / 7, 0.3, a, returnt = 1)
Mct_a10_03m = ph.t_vs_Mc_fixed_a(2. / 5, 0.3, a, returnt = 1)

Mc00 = Mct_a10_00[0]
ttot00 = Mct_a10_00[1]
Mcrittot00 = Mct_a10_00[2]
dt00 = Mct_a10_00[3]
tcum00 = Mct_a10_00[4]

Mc03 = Mct_a10_03[0]
ttot03 = Mct_a10_03[1]
Mcrittot03 = Mct_a10_03[2]
dt03 = Mct_a10_03[3]
tcum03 = Mct_a10_03[4]

Mc03m = Mct_a10_03m[0]
ttot03m = Mct_a10_03m[1]
Mcrittot03m = Mct_a10_03m[2]
dt03m = Mct_a10_03m[3]
tcum03m = Mct_a10_03m[4]


savefig = 0
filename = 'coolingtime_vs_Mc_10au.pdf'

plt.semilogy(Mc00, ttot00 / 10**6, '-o', label = r'$\nabla_{\mathrm{ad}}=2/7$, $\mu=$' + str(mufn(Y00)))
plt.semilogy(Mc03, ttot03 / 10**6, '-sr', label = r'$\nabla_{\mathrm{ad}}=2/7$,  $\mu=$' + str(mufn(Y03))[:4])
plt.semilogy(Mc03m, ttot03m / 10**6, '-g^', label = r'$\nabla_{\mathrm{ad}}=2/5$,  $\mu=$' + str(mufn(Y03))[:4])
plt.axhline(y = 3, linestyle = '--', color = 'black', label = 't = 3 Myrs')
plt.xlabel(r'$M_c (M_{\bigoplus})$ ')
plt.ylabel(r'$t_{\mathrm{cool}}$(Myrs)')
#plt.xlim(0, 15)
plt.legend()
plt.title(r'a=' + str(int(a)) + ' AU')
if savefig != 0:
    plt.savefig(userpath + '/figs/ModelAtmospheres/RadSelfGravPoly/' + filename)
plt.show()

au = [1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
#au = numpy.linspace(10, 100, 10)

crit00 = ph.Mc_vs_a_fixed_t(2. / 7, 0.0, 3 * 10**6)
crit03 = ph.Mc_vs_a_fixed_t(2. / 7, 0.3, 3 * 10**6)
crit03m = ph.Mc_vs_a_fixed_t(2. / 5, 0.3, 3 * 10**6)

au00 = crit00[0]
au03 = crit03[0]
au03m = crit03m[0]

Mcrit00 = crit00[1]
Mcrit03 = crit03[1]
Mcrit03m = crit03m[1]



fig = plt.figure() # (figsize = (14.17/2, 6.04/2))
ax = fig.add_subplot(1,1,1)

savefig = 1

ax.semilogx(au00, Mcrit00, '-o', label = r'$\nabla_{\mathrm{ad}}=2/7$, $\mu=$' + str(mufn(Y00))[:4])
ax.semilogx(au03, Mcrit03, '-rs', label = r'$\nabla_{\mathrm{ad}}=2/7$, $\mu=$' + str(mufn(Y03))[:4])
ax.semilogx(au03m, Mcrit03m, '-g^', label = r'$\nabla_{\mathrm{ad}}=2/5$, $\mu=$' + str(mufn(Y03))[:4])
ax.set_xlabel('a [AU]')
ax.set_ylabel(r'$M_c [M_{\oplus}]$ ')
ax.set_xticks((1, 5, 10, 20, 40, 60, 80, 100))
ax.set_xticklabels(("1", "5", "10", "20", "40", "60", "80", "100"))
ax.set_xlim(0.8, 110)
#ax.set_ylim(3, 12)
ax.legend()
ax.set_title(r'Critical core mass vs. distance for disk lifetime of 3 Myrs')
if savefig == 1:
    plt.savefig(userpath + \
                '/figs/ModelAtmospheres/RadSelfGravPoly/Mcrit_vs_a_3Myrs_new3.pdf')
plt.show()




