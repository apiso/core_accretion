import os
os.chdir('/home/apiso/aacore/python')

from RadSGPoly.atmseries_poly import atmseries
from utils.constants import Me

a = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]

for i in range(len(a)):
    atmseries(a[i], 3.2, 2./7, 0.3, 1.0*Me, 10.0*Me, 10)
