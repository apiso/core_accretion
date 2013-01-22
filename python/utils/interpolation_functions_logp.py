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



EOS=EOSread.EOSread(Nt, Np, Ncol, userpath + '/dat/EOStables/H_HE_TAB_I_30_EXT.txt')[0]
T=EOSread.EOSread(Nt, Np, Ncol, userpath + '/dat/EOStables/H_HE_TAB_I_30_EXT.txt')[1]
EOSs=EOSread.EOSread(Np, 25, 6, userpath + '/dat/EOStables/EOS_S_H_HE_30_PRESSURE.txt')[0]

p = EOS[-1, :, 0]
x = p

#entropy interpolation from temperature and pressure
y = EOS[-1, :, 0]
z = EOS[:, :, 4]
yy, xx = numpy.meshgrid(y, T)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolations = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(t), log(p)


#delad interpolation from pressure and entropy
y = EOSs[-1,:,0]
z = EOSs[:,:,3]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolationdelad = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(p), s


#fact (= P/(rho*delad)) interpolation from pressure and entropy
y = EOSs[-1,:,0]
z = EOSs[:,:,2]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolationfact = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(p), s


#temperature interpolation from pressure and entropy
y = T
z = EOS[:, :, 4]
xx, yy = numpy.meshgrid(x, y)
pts = numpy.array((xx.ravel(), z.ravel())).T
sr = yy.ravel()
interpolationt = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(p), s


#density interpolation from pressure and entropy
y = EOSs[-1,:,0]
z = EOSs[:,:,4]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolationrho = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(p), s


#internal energy interpolation from pressure and density
y = EOSs[-1,:,0]
z = EOSs[:,:,5]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolationu = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(p), s


