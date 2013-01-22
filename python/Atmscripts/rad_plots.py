import numpy
from RadNoSGPoly import write_read_poly as wrp
from RadNoSGRealGas import write_read_realgas as wrr
from RadSGPoly import wr_selfgrav_poly as wrsgp
from RadNoSGPoly import atmprofiles_poly as atm
from RadSGPoly import atmprofiles_selfgrav as atmsg
from utils import userpath
import matplotlib.pyplot as plt

##profpolylow = wrp.atmprofileread(userpath.userpath + '/dat/MODELS/RadNoSGPoly/Mc5_low_prof.txt')[0]
##profpolyhigh = wrp.atmprofileread(userpath.userpath + '/dat/MODELS/RadNoSGPoly/Mc5_high_prof.txt')[0]
##parampolylow = wrp.atmparamread(userpath.userpath + '/dat/MODELS/RadNoSGPoly/Mc5_low_param.txt')[0]
##parampolyhigh = wrp.atmparamread(userpath.userpath + '/dat/MODELS/RadNoSGPoly/Mc5_high_param.txt')[0]
##Spoly = wrp.atmparamread(userpath.userpath + '/dat/MODELS/RadNoSGPoly/Mc5_high_param.txt')[1]
##
##profreallow = wrr.atmprofileread(userpath.userpath + '/dat/MODELS/RadNoSGRealGas/Mc5_low_prof.txt')[0]
##profrealhigh = wrr.atmprofileread(userpath.userpath + '/dat/MODELS/RadNoSGRealGas/Mc5_high_prof.txt')[0]
##paramreallow = wrr.atmparamread(userpath.userpath + '/dat/MODELS/RadNoSGRealGas/Mc5_low_param.txt')[0]
##paramrealhigh = wrr.atmparamread(userpath.userpath + '/dat/MODELS/RadNoSGRealGas/Mc5_high_param.txt')[0]
##Sreal = wrr.atmparamread(userpath.userpath + '/dat/MODELS/RadNoSGRealGas/Mc5_high_param.txt')[1]
##
##profpolysg = wrsgp.atmprofileread(userpath.userpath + '/dat/MODELS/RadSGPoly/Mc5_prof.txt')
##parampolysg = wrsgp.atmparamread(userpath.userpath + '/dat/MODELS/RadSGPoly/Mc5_param.txt')
##
##Lpolylow = parampolylow[:,-1]
##Lpolyhigh = parampolyhigh[:,-1]
##Lreallow = paramreallow[:,-1]
##Lrealhigh = paramrealhigh[:,-1]
##Lpolysg = parampolysg[:,-3]
##
##Mpolylow = parampolylow[:,0]
##Mpolyhigh = parampolyhigh[:,0]
##Mradpolylow = parampolylow[:,1]
##Mradpolyhigh = parampolyhigh[:,1]
##Mradreallow = paramreallow[:,1]
##Mradrealhigh = paramrealhigh[:,1]
##Mreallow = paramreallow[:,0]
##Mrealhigh = paramrealhigh[:,0]
##Mpolysg = parampolysg[:,0]
##Mradpolysg = Mpolysg - parampolysg[:,1]
##
##Mpoly = numpy.append(Mpolylow[::-1], Mpolyhigh)
##Mradpoly = numpy.append(Mradpolylow[::-1], Mradpolyhigh)
##Lpoly = numpy.append(Lpolylow[::-1], Lpolyhigh)
##Spoly = numpy.append(Spoly[::-1], Spoly)
##
##Mreal = numpy.append(Mreallow[::-1], Mrealhigh)
##Mradreal = numpy.append(Mradreallow[::-1], Mradrealhigh)
##Lreal = numpy.append(Lreallow[::-1], Lrealhigh)
##Sreal = numpy.append(Sreal[::-1], Sreal)

plt.figure(1)
plt.semilogx(Mpoly, Spoly, 'b', linewidth = 2)
plt.semilogx(Mreal, Sreal, '--r', linewidth = 2)
plt.legend(('Polytrope', 'Real EOS'))
plt.xlabel(r'$M$ ($M_{\oplus}$)')
plt.ylabel('log (S)')
plt.title(r'$a=10$ AU, $M_c=5 M_{\oplus}$')
plt.savefig(userpath.userpath + '/figs/ModelAtmospheres/Rad_Compare/smplotMc5_2.pdf')
plt.show()

##f, sm = plt.subplots(2, sharex = True)
##sm[0].semilogx(Mpolylow, Spoly, Mpolyhigh, Spoly, 'b', linewidth = 2)
##sm[0].set_ylabel('log (S)')
##sm[0].set_title(r'Polytrope no SG: $a=10$ AU, $M_c=5 M_{\oplus}$')
##sm[1].semilogx(Mreallow, Sreal, Mrealhigh, Sreal, 'b', linewidth = 2)
##sm[1].set_xlabel(r'$M$ ($M_{\oplus}$)')
##sm[1].set_ylabel('log (S)')
##sm[1].set_title(r'Real EOS no SG: $a=10$ AU, $M_c=5 M_{\oplus}$')
##plt.savefig(userpath.userpath + '/figs/ModelAtmospheres/Rad/smplotMc5.pdf')
##plt.show()

plt.figure(2)
plt.loglog(Mpoly, Lpoly, 'b', linewidth = 2)
plt.loglog(Mreal, Lreal, '--r', linewidth = 2)
plt.loglog(Mpolysg, Lpolysg, '-.g', linewidth = 2)
plt.xlabel(r'$M$ ($M_{\oplus}$)')
plt.ylabel(r'$L$ (erg s$^{-1}$)')
plt.title(r'$a=10$ AU, $M_c=5 M_{\oplus}$')
plt.legend(('Polytrope no SG', 'Real EOS no SG', 'Polytrope SG'), loc = 2)
plt.savefig(userpath.userpath + '/figs/ModelAtmospheres/Rad_Compare/lmplotMc5_2.pdf')
plt.show()

plt.figure(3)
plt.loglog(Mpoly, Mradpoly, 'b', linewidth = 2)
plt.loglog(Mreal, Mradreal, '--r', linewidth = 2)
plt.loglog(Mpolysg, Mradpolysg, '-.g', linewidth = 2)
plt.xlabel(r'$M$ ($M_{\oplus}$)')
plt.ylabel(r'$M_{\mathrm{rad}}$ ($M_{\oplus}$)')
plt.title(r'Real EOS no SG: $a=10$ AU, $M_c=5 M_{\oplus}$')
plt.legend(('Polytrope no SG', 'Real EOS no SG', 'Polytrope SG'), loc = 2)
plt.savefig(userpath.userpath + '/figs/ModelAtmospheres/Rad_Compare/mradmplotMc5_2.pdf')
plt.show()

##f2, lm = plt.subplots(3, sharex = True)
##lm[0].semilogy(Mpolylow, Lpolylow, Mpolyhigh, Lpolyhigh, 'b', linewidth = 2)
##lm[0].set_ylabel(r'$L$ (erg s$^{-1}$)')
##lm[0].set_title(r'Polytrope no SG: $a=10$ AU, $M_c=5 M_{\oplus}$')
##lm[1].semilogy(Mreallow, Lreallow, Mrealhigh, Lrealhigh, 'b', linewidth = 2)
##lm[1].set_ylabel(r'$L$ (erg s$^{-1}$)')
##lm[1].set_title(r'Real EOS no SG: $a=10$ AU, $M_c=5 M_{\oplus}$')
##lm[2].semilogy(Mpolysg, Lpolysg, 'b', linewidth = 2)
##lm[2].set_ylabel(r'$L$ (erg s$^{-1}$)')
##lm[2].set_title(r'Polytrope SG: $a=10$ AU, $M_c=5 M_{\oplus}$')
##lm[2].set_xlabel(r'$M$ ($M_{\oplus}$)')
##lm[2].set_xlim(0,100)
##plt.savefig(userpath.userpath + '/figs/ModelAtmospheres/Rad/lmplotMc5.pdf')
##plt.show()
##
##f3, mm = plt.subplots(3, sharex = True)
##mm[0].plot(Mpolylow, Mradpolylow, Mpolyhigh, Mradpolyhigh, 'b', linewidth = 2)
##mm[0].set_ylabel(r'$M_{\mathrm{rad}}$ ($M_{\oplus}$)')
##mm[0].set_title(r'Polytrope no SG: $a=10$ AU, $M_c=5 M_{\oplus}$')
##mm[1].plot(Mreallow, Mradreallow, Mrealhigh, Mradrealhigh, 'b', linewidth = 2)
##mm[1].set_ylabel(r'$M_{\mathrm{rad}}$ ($M_{\oplus}$)')
##mm[1].set_title(r'Real EOS no SG: $a=10$ AU, $M_c=5 M_{\oplus}$')
##mm[2].plot(Mpolysg, Mradpolysg, 'b', linewidth = 2)
##mm[2].set_ylabel(r'$M_{\mathrm{rad}}$ ($M_{\oplus}$)')
##mm[2].set_title(r'Polytrope SG: $a=10$ AU, $M_c=5 M_{\oplus}$')
##mm[2].set_xlabel(r'$M$ ($M_{\oplus}$)')
##plt.savefig(userpath.userpath + '/figs/ModelAtmospheres/Rad/mradmplotMc5.pdf')
##plt.show()


plt.subplot(2,1,1)
atm.mrplot(profpolylow[0:-1:5])
atm.mrplot(profpolyhigh[0:-1:5])
plt.title(r'Polytrope no SG: $a=10$ AU, $M_c=5 M_{\oplus}$', fontsize = 10)
#plt.text(10**(-4), 2*10**4, 'Polytrope no SG')
plt.subplot(2,1,2)
atm.mrplot(profreallow[0:-1:3])
atm.mrplot(profrealhigh[0:-1:3])
#plt.text(10**(-4), 5*10**4, 'Real EOS no SG')
plt.xlabel(r'$M_{\mathrm{atm}}$ ($M_{\oplus}$)')
plt.title(r'Real EOS no SG: $a=10$ AU, $M_c=5 M_{\oplus}$', fontsize = 10)
#plt.savefig(userpath.userpath + '/figs/ModelAtmospheres/Rad/mrprofileplotMc5.pdf')
plt.show()







