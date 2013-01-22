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
EOSs=EOSread.EOSread(Nt, 25, 6, userpath + '/dat/EOStables/EOS_S_H_HE_30.txt')[0]


x = T

#entropy interpolation from temperature and pressure
y = EOS[-1, :, 0]
z = EOS[:, :, 4]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolations = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(t), log(p)


#delad interpolation from temperature and entropy
y = EOSs[-1,:,0]
z = EOSs[:,:,3]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolationdelad = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(t), s


#fact (= P/(rho*delad)) interpolation from temperature and entropy
y = EOSs[-1,:,0]
z = EOSs[:,:,2]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolationfact = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(t), s


#pressure interpolation from temperature and entropy
y = EOSs[-1,:,0]
z = EOSs[:, :, 1]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolationp2 = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(t), s


#density interpolation from temperature and entropy
y = EOSs[-1,:,0]
z = EOSs[:,:,4]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolationrho = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(t), s


#internal energy interpolation from temperature and density
y = EOSs[-1,:,0]
z = EOSs[:,:,5]
yy, xx = numpy.meshgrid(y, x)
pts = numpy.array((xx.ravel(), yy.ravel())).T
sr = z.ravel()
interpolationu = scipy.interpolate.LinearNDInterpolator(pts,sr) #args: log(t), s


