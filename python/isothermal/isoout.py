import numpy as np
from isotherm import G, R , Po, To , Ho , Me , a #since this module is in same folder there should be no problem finding...
from isotherm import isostructlnr, odeint , Mdisk , RHill
from scipy.interpolate import interp1d

mc = 5*Me
fc = 0.1 #location of base of isothermal layer relative to core's Bondi radius.

def rBc(mc, To = To):
    return G * mc / R / To

def intout(Pover, mc = mc,  sg = 1, ngrid = 200 , retall = 0):
    """
    outward integration with a base (core) pressure that is "Pover" times the disk pressure.
    This fn. returns the radius and mass where the pressure drops to the disk pressure 
    """
    rBc = G * mc / R / To
    rc = fc * rBc  #puts "core",  i.e. base of isothermal layer, at some fraction of core's Bondi radius.
    lnrgrid = np.linspace(np.log(rc) , np.log(Ho) , ngrid)
    f , bc = isostructlnr , [mc, np.log(Pover * Po) , np.log(To)]
    y = odeint(f, bc, lnrgrid, args = (sg,mc)) #integration of the structure equations
    p = np.exp(y[ : , 1])
    if p[-1] > Po:
        print "atmos. doesn't drop to disk pressure at log10(Pover) = ", np.log10(Pover)
        import sys
        sys.exit()
    lnp_rev = y[::-1 , 1] #log pressure in increasing order
    lnr_rev = lnrgrid[::-1] #must reverse other arrays to match (stupid interp1d...)
    m_rev = y[::-1 , 0]
    f = interp1d(lnp_rev , lnr_rev)
    rmatch = np.exp(f(np.log(Po)))
    fmp = interp1d(lnp_rev , m_rev)
    mmatch = fmp(np.log(Po))
    if retall:
        return rmatch/rBc , mmatch/mc , y , np.exp(lnrgrid)#y& r for testing
    else:
        return rmatch/rBc , mmatch/mc 

def varyPbase(LP1,LP2,nPr, mc = mc, sg = 1):
    """
    parameters:
    LP1 , LP2 -- log10 of range of base pressures (e.g.) 2 is an overdensity of 100.
    nPr --- number of solutions (log spacing in pressure used)

    returns:
    ra --- matching radius in terms of CORE bondi radius
    ma --- matching mass (including disk mass) in terms of core mass
    
    """
    Povera = np.logspace(LP1, LP2, nPr)
    ra , ma = np.zeros(nPr) , np.zeros(nPr)
    for i , Pover in enumerate(Povera):
        (ra[i] , ma[i]) = intout(Pover, mc , sg)
    return ra, ma , Povera

def mnorms(ma , ra , mc = mc):
    """disk correction to atmosphere's mass"""
    rBc = G * mc / R / To ; rc = fc * rBc
    Mdiska = Mdisk(ra * rBc , rc)
    mcorr = ma - Mdiska/mc
    return mcorr

def rnorms(ma, ra , mc = mc):
    """
    input: ra --- rmatch over core's bondi radius
    returns: (rm_Bt , rm_Bds , rm_Hc , rm_Ht , rm_Hds)
    arrays of rmatch normalized to total (Bt) and total disk corrected (Bds) Bondi radius,
       and Hill radius for core (Hc),  total (Ht) and disk sbtracted (Hds) masses.
    """
    mtot = ma * mc
    mcorr = mnorms(ma , ra , mc) * mc # in cgs
    rBc = G * mc / R / To
    rBt = rBc * mtot / mc
    rBds = rBc * mcorr / mc
    rHc = RHill(mc , a)
    rHt , rHds = (mtot / mc)**(1./3) * rHc , (mcorr / mc)**(1./3) * rHc 
    rout = ra * rBc
    return (rout / rBt , rout / rBds , rout / rHc  , rout / rHt  , rout / rHds)
