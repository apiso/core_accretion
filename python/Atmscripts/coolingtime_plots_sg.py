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

delad = 2./5


load = 1
cool = 1
savefig = 0
nofixeddisk = 0

if nofixeddisk == 0:
    Y = 0.3
    a = 1.0

prms = params(Me, Re, a, delad, Y, gamma = gammafn(delad), R = Rfn(Y), \
              Cv = Cvfn(Y, delad), Pd = Pdisk(a, mstar, FSigma, FT), \
              Td = Tdisk(a, FT), kappa = kdust)

if load != 0:

    if Y == 0.3:
	if delad == 2./7 and a == 5.0:
		atm1 = atmload('Mc7.0.npz', prms = prms)
		atm2 = atmload('Mc8.0.npz', prms = prms)
		atm3 = atmload('Mc9.0.npz', prms = prms)
		atm4 = atmload('Mc10.0.npz', prms = prms)
		atm5 = atmload('Mc11.0.npz', prms = prms)

		Mc = numpy.linspace(7, 11, 5)

	elif delad == 2./7 and a == 1.0:
		atm1 = atmload('Mc10.0.npz', prms = prms)
		atm2 = atmload('Mc11.0.npz', prms = prms)
		atm3 = atmload('Mc12.0.npz', prms = prms)
		atm4 = atmload('Mc13.0.npz', prms = prms)
		atm5 = atmload('Mc14.0.npz', prms = prms)

		Mc = numpy.linspace(10, 14, 5)
	
	elif delad == 2./5 and a == 5.0:
		atm1 = atmload('Mc6.0.npz', prms = prms)
		atm2 = atmload('Mc7.0.npz', prms = prms)
		atm3 = atmload('Mc8.0.npz', prms = prms)
		atm4 = atmload('Mc9.0.npz', prms = prms)
		atm5 = atmload('Mc10.0.npz', prms = prms)

		Mc = numpy.linspace(6, 10, 5)

	elif delad == 2./5 and a == 1.0:
		atm1 = atmload('Mc30.0.npz', prms = prms)
		atm2 = atmload('Mc39.0.npz', prms = prms)
		atm3 = atmload('Mc39.2.npz', prms = prms)
		atm4 = atmload('Mc39.5.npz', prms = prms)
		atm5 = atmload('Mc45.0.npz', prms = prms)

		Mc = [26.0, 27.0, 39.0, 39.25, 39.5] #numpy.linspace(16, 20, 5)

	else:


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

        	Mc = numpy.linspace(1, 10, 10)

    elif Y == 0.0:

	if delad == 2./7 and a == 5.0:
		atm1 = atmload('Mc11.0.npz', prms = prms)
		atm2 = atmload('Mc12.0.npz', prms = prms)
		atm3 = atmload('Mc13.0.npz', prms = prms)
		atm4 = atmload('Mc14.0.npz', prms = prms)
		atm5 = atmload('Mc15.0.npz', prms = prms)

		Mc = numpy.linspace(11, 15, 5)

	elif delad == 2./7 and a == 1.0:
		atm1 = atmload('Mc16.0.npz', prms = prms)
		atm2 = atmload('Mc17.0.npz', prms = prms)
		atm3 = atmload('Mc18.0.npz', prms = prms)
		atm4 = atmload('Mc19.0.npz', prms = prms)
		atm5 = atmload('Mc20.0.npz', prms = prms)

		Mc = numpy.linspace(16, 20, 5)

	else:
        
        	atm5 = atmload('Mc5.0.npz', prms = prms)
        	atm6 = atmload('Mc6.0.npz', prms = prms)
        	atm7 = atmload('Mc7.0.npz', prms = prms)
        	atm8 = atmload('Mc8.0.npz', prms = prms)
        	atm9 = atmload('Mc9.0.npz', prms = prms)
        	atm10 = atmload('Mc10.0.npz', prms = prms)
        	atm11 = atmload('Mc11.0.npz', prms = prms)
        	atm12 = atmload('Mc12.0.npz', prms = prms)
        	atm13 = atmload('Mc13.0.npz', prms = prms)
        	atm14 = atmload('Mc14.0.npz', prms = prms)       

        	Mc = numpy.linspace(5, 14, 10)

if a == 1.0 or a == 5.0:
	model = [atm1[0], atm2[0], atm3[0], atm4[0], atm5[0]]
	param = [atm1[1], atm2[1], atm3[1], atm4[1], atm5[1]]
	prof = [atm1[2], atm2[2], atm3[2], atm4[2], atm5[2]]

else:
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
plt.title(r'$\nabla=$' + str(delad) + ', Y=' + str(Y)[:4] + ', a=' + str(int(a)) + ' AU')
if savefig != 0:
    plt.savefig(userpath + '/figs/ModelAtmospheres/RadSelfGravPoly/coolingtime_vs_Mc_' + str(int(a)) + 'au.pdf')
plt.show()



