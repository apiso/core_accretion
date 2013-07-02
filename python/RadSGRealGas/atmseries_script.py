import os
#os.chdir('/home/apiso/repos/core_accretion/python')
os.chdir('/Users/ana-mariapiso/core_accretion/python')

import numpy
from RadSGRealGas.atmseries import atmseries
from RadSGRealGas.profiles import atmload
from RadSGRealGas.cooling import critical as crit
from utils.constants import Me, kfixed, kdustbeta1, \
     Pdisk, Tdisk, kdust, kdustall, paramsEOS
from utils.parameters import mstar, FSigma, FT, rhoc

a = [100.0] # numpy.append([5.0], numpy.linspace(10, 100, 10))
MBi = 0 * numpy.ndarray(shape = (len(a)), dtype = float)
MBf = 0 * numpy.ndarray(shape = (len(a)), dtype = float)
t = 0 * numpy.ndarray(shape = (len(a)), dtype = float)


for i in range(len(a)):
    Mc = 3 * Me
    rc = (3 * Mc / (4 * numpy.pi * rhoc))**(1./3)
    prms = paramsEOS(Mc, rc, 0.3, a[i], Pd=Pdisk(a[i], mstar, FSigma, FT), \
                   Td=Tdisk(a[i],FT), kappa=kfixed)
    x = atmload('Mc3.0.npz', prms=prms)
    model, param, prof = x
    y = crit(param, prof, model)

    MBi[i] = y[0].MB[0]
    MBf[i] = y[0].MB[-1]
    t[i] = sum(y[2])/(365*24*3600*10**6)
    
    

##for i in range(len(a)):
##    atmseries(a[i], 3.2, 0.3, 3.0*Me, 5.0*Me, 3, nMpoints = 100,  minMfrac = 1.05, \
##              maxMfrac = 2.3, opacity = kdustbeta1)
####
#for i in range(len(a)):
#    atmseries(a[i], 3.2, 2./7, 0.0, 16.0*Me, 20.0*Me, 5, minMfrac = 1.2, maxMfrac = 2.5)

#for i in range(len(a)):
#    atmseries(a[i], 3.2, 2./5, 0.3, 15.0*Me, 19.0*Me, 5, minMfrac = 1.05, maxMfrac = 2.5)
