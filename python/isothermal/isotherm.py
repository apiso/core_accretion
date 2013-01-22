"""
module for integrating structure of isothermal atmosphere that matches onto disk.  Goal is to check for convergence with different choices of nfit.  The disk mass is no longer included as part of the atmosphere mass  

usage (from python folder or with userpath loaded)
>>> from isothermal import isotherm as iso #(abbreviation optional)
>>> iso.trialint(mguess, nfit, sg)  # gives fractional deviation of core mass for a trial mass (lnr integration)
>>> iso.trialintm(mguess, nfit) # gives deviation of core radius (m integration) 
>>> iso.trymanym(params...) #get fractional fit for a range of mass guesses to look for a root
>>> iso.P #pressure vector of trial solution (and other quantities)
"""

import numpy as np
from constants import G , Me , Tdisk , Pdisk , Hdisk , RHill
from parameters import R, rhoc
from scipy.integrate import odeint
from scipy.optimize import bisect
import sys

#model parameters (could import more of these)
a = 10;
mc = 5* Me
Po , To = Pdisk(a , 1. , 1. , 1.) , Tdisk (a, 1.)
Ho = Hdisk(a , 1.)

#derived model params
rhoo = Po/(R * To)
#rc =  (3*mc/(4*np.pi*rhoc))**(1./3)
rBc = G * mc / (R * To)
rc = 0.1 * rBc  #puts "core" at some fraction of core's Bondi radius, to mimic reality that isothermal layer does not extend to physical core

def rfit_fn(nfit , mg = mc, a = a):
    """
    gives fit radius as Hill radius (for nfit = 0) or some multiple (or fraction) nfit <> 0 of Bondi radii
    """
    if (nfit == 0) or (nfit == 'Hill'):
        rfit =  RHill(mg , a)
    else:
        rfit =  nfit * G * mg / (R * To)
    return min(rfit , Ho)

def Mdisk(rout, rc = rc):
    """
    returns mass of uncompressed disk gas in "atmosphere"
    """
    return (4 * np.pi / 3) * rhoo * (rout**3 - rc**3)

def isostructm(y, m):
    """
    structure eqns. dy/dm for "y" = [ln(r) , ln(P), ln(T)]
    dln(r)/dm = 1 / (4 * pi * r**3 * rho),
    dln(P)/dm = - G * m / (4 * pi * r**4 * P)
    dln(T)/dm = 0
  
    """
    rho = np.exp(y[1] - y[2]) / R  
    return np.array([ (4 * np.pi * rho * np.exp(3 * y[0]))**(-1), \
                         - G * m / (4 * np.pi * np.exp(4 * y[0] + y[1])), \
                         0])

def isostructlnr(y, lnr, sg = 1, mco = mc):
    """
    structure eqns. dy/dlnr for "y" = [m, ln(P), ln(T)]
    dm/dln(r) = (4 * pi * r**3 * rho),
    dln(P)/dln(r) = - G * m * rho / r**3
    dln(T)/dln(r) = 0

    Keywords:
    sg -- include self-gravity in hydrostatic balance (default 1 = True)
    """
    if sg:
        m = y[0]
    else:
        m = mco
    rho = np.exp(y[1] - y[2]) / R   
    return np.array([ 4 * np.pi * np.exp(3 * lnr) * rho, \
                         - G * m *rho / np.exp(lnr + y[1]), \
                         0])

def isostructlnr_nsg(y, lnr):
    """
    structure eqns. dy/dlnr for "y" = [m, ln(P), ln(T)]
    dm/dln(r) = (4 * pi * r**3 * rho),
    dln(P)/dln(r) = - G * m_core * rho / r**3 
    dln(T)/dln(r) = 0
    """
    rho = np.exp(y[1] - y[2]) / R   
    return np.array([ 4 * np.pi * np.exp(3 * lnr) * rho, \
                         - G * mc *rho / np.exp(lnr + y[1]), \
                         0])

def trialintm(mguess, nfit = 0):
    """
    integrate structure equations vs. mass from top and compare core radius of solution to actual value.
    Returns log of core radius ratio, which crosses zero at the root.
    Note that mguess does NOT include the disk mass correction, so it is the physical mass we are interested in.
    """
    rfit = rfit_fn(nfit , mguess)
    f , bc = isostructm , [np.log(rfit), np.log(Po) , np.log(To)]
    mout = mguess + Mdisk(rfit)
    y = odeint(f, bc, [mout , mc]) #integration of the structure equations
    return y[-1 , 0] - np.log(rc)

def trialint(mguess, nfit, sg = 0):
    """
    Integrate structure equations vs lnr from top and compare core mass of solution to actual value. 
    Returns core mass ratio - 1, which crosses zero at the root.
    Note that mguess does NOT include the disk mass correction, so it is the physical mass we are interested in.
    """
    #rfit = nfit * G * mc / (R * To) #possible problems with using mguess...
    rfit = rfit_fn(nfit , mguess)
    mout = mguess + 0*Mdisk(rfit)
    f , bc = isostructlnr , [mout, np.log(Po) , np.log(To)]    
    y = odeint(f, bc, [np.log(rfit) , np.log(rc)], args = (sg,)) #integration of the structure equations
    return y[-1 , 0]/mc - 1    

def trymanym(mmax, mmin = mc, sg = 1, nm = 20, nfit = 4,plotit = 0):
    """
    return closeness of fit for a range of masses (evenly spaced)
    optional plotting option
    """
    ms , err = np.linspace(mmin,mmax,nm) , np.zeros(nm)
    for i, mguess in enumerate(ms):
        err[i] = trialint(mguess , nfit, sg)
    if plotit:
        import matplotlib.pyplot as plt
        plt.plot(ms/Me , err)
    return err

def zbrak(nfit, sg = 1, mc = mc, ma1 = .1*mc , ma2 = 0.2*mc , Po = Po, To = To):
    """
    attempt to braket solution for atmospheric mass, returning braketed values
    #TODO add testing for local min/max
    """
    ntry , factor = 50 , 1.6
    #print sg
    err1 = trialint(mc + ma1, nfit, sg = sg)
    err2 = trialint(mc + ma2, nfit, sg = sg)
    for t in range(ntry):
        #print err1,err2
        if (err1 * err2 < 0):
            print "atm mass between", ma1 / Me, ma2 / Me
            return (ma1 , ma2 ) , (err1 , err2)
        else:
            if (abs(err1) < abs(err2)):
                ma1 = ma1 + factor * ( ma1 - ma2 )
                err1 = trialint(mc + ma1, nfit, sg = sg)
            else:
                ma2 = ma2 + factor * ( ma2 - ma1 )
                err2 = trialint(mc + ma2, nfit, sg = sg)

    print "Failure to bracket roots between masses (in Earth masses)", ma1/Me, ma2/Me
    sys.exit()       

def atmmass(nfit, sg = 1, mc = mc, ma1 = .1*mc , ma2 = 0.2*mc , Po = Po, To = To):    
    """
    loop to find atmosphere mass (in progress)
    """
    maxtries , tol , t  = 50, 1e-5 , 1
    #bracketed solution
    (ma1 , ma2) , (err1 , err2) = zbrak(nfit, sg, mc = mc, ma1 = .1*mc , ma2 = 0.2*mc)
    errmin = min(abs(err1) , abs(err2))
    ms = bisect(trialint , mc + ma1 , mc + ma2 , args = (nfit,sg) )

    print "root found with atm mass", (ms - mc)/Me, "and error", trialint(ms , nfit , sg)
    return ms - mc
    
    # while ( (errmin > tol) and (try < maxtries) ):
    #     try +=1
    #     if newton == 1:
    #         #take newton step
            
  
    return (err1, err2)



#sample integrations for command line evaluation
ngrid , nfit = 100 , 0.
mi = 6.6*Me #initial guess for isothermal atmosphere + core mass
rfit = rfit_fn(nfit, mi)

#radius grid
lnrgrid = np.linspace(np.log(rfit) , np.log(rc) , ngrid)
r = np.exp(lnrgrid)
mout = mi + 0*Mdisk(rfit)
f , bc = isostructlnr , [mout, np.log(Po) , np.log(To)]    
y = odeint(f, bc, lnrgrid)
m , P = y[: , 0]/Me , np.exp(y[:,1])
y_nsg = odeint(f, bc, lnrgrid, args = (0,))
m_nsg , P_nsg = y_nsg[: , 0]/Me , np.exp(y_nsg[:,1])
P_an = Po * np.exp( rBc / np.exp(lnrgrid) - rBc / rfit )

#mass grid
mgrid = np.linspace(mi + Mdisk(rfit) , mc , ngrid) #descending in mass
f , bc = isostructm , [np.log(rfit), np.log(Po) , np.log(To)]
ym = odeint(f, bc, mgrid) #integration of the structure equations
rm , Pm = np.exp(ym[: , 0]) , np.exp(ym[:,1])
rcint = np.exp(ym[-1 , 0])

#parameter for evaluating uncompressed disk mass
theta = 4 * np.pi * rhoo * rBc**3 / (3 * mc)
