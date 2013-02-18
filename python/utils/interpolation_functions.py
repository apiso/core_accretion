"""

This module defines the 2-D interpolation functions which will then be applied to pairs of
thermodynamic variables; used in the shooting code, shooting_odeint

"""

import numpy
import scipy
from scipy import interpolate
from EOS import EOSread
from table_size import Nt, Np, Ncol
from userpath import userpath



EOS=EOSread.EOSread(Nt, Np, Ncol, userpath + \
                    '/dat/EOStables/H_HE_TAB_I_30_EXT.txt')[0]
log10T=EOSread.EOSread(Nt, Np, Ncol, userpath + \
                  '/dat/EOStables/H_HE_TAB_I_30_EXT.txt')[1]

log10p = EOS[-1, :, 0] #log pressure

log10rho = EOS[:, :, 3] #log density
log10delad = EOS[:, :, 10]
log10u = EOS[:, :, 5]
log10S = EOS[:, :, 4]

#density interpolation from temperature and pressure
yy, xx = numpy.meshgrid(log10p, log10T)
pts = numpy.array((xx.ravel(), yy.ravel())).T

srrho = log10rho.ravel()
srd = log10delad.ravel()
sru = log10u.ravel()
srs = log10S.ravel()

interplog10rho = scipy.interpolate.LinearNDInterpolator(pts,srrho)
interpdelad = scipy.interpolate.LinearNDInterpolator(pts,srd)
interplog10u = scipy.interpolate.LinearNDInterpolator(pts,sru)
interplog10S = scipy.interpolate.LinearNDInterpolator(pts,srs) \
                   #args: log(t), log(p)

