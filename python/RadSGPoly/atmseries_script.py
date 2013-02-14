import os
os.chdir('/home/apiso/repos/core_accretion/python')


from RadSGPoly.atmseries_poly import atmseries
from utils.constants import Me

a = [1.75]

#for i in range(len(a)):
#    atmseries(a[i], 3.2, 2./7, 0.3, 12*Me, 14.0*Me, 3, minMfrac = 1.2, maxMfrac = 2.5)

#for i in range(len(a)):
#    atmseries(a[i], 3.2, 2./7, 0.0, 16.0*Me, 20.0*Me, 5, minMfrac = 1.2, maxMfrac = 2.5)

for i in range(len(a)):
    atmseries(a[i], 3.2, 2./5, 0.3, 15.0*Me, 19.0*Me, 5, minMfrac = 1.05, maxMfrac = 2.5)
