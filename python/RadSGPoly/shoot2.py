"""
Find atmosphere solutions for a polytrope EOS and fully self-gravitating
atmosphere out to Hill radius.  Both `static` (fixed luminosity) and `evolving`
solutions.

Loosely there are (at least) 6 levels of founctions here:
0) physics helpers: `delRad`, interpolation fns (generated by `makeprev`)
1) structure eqns.: `polyradstruct_r` (static or evo, but no conv. switch)
                    `evostruct_r` (evo only, conv. switch)
2) integrators: we use scipy.optimize.odeint
3) odeint wrappers: integrate structure eqns. from top to bottom for:
    a) Error calculation: `error1D`, `error2D`
    b) Generating a profile or relevant info of a converged solution:
       `profiles` (1D), `makeprev` (1D or 2D), staticguess (1D but used as guess
       for 2D rootfinders.  (One could argue that 3b should be a higher level
       than 4 b/c rootfinding results are used.  Whatever.)
4) rootfinders: `lLtop` (1D) `root2D`.  Functions which map the error space
                without rootfinding are `errorrange` (1D) and `errormap` (2D).
                These are simple wrappers of scipy.optimize rootfinders.
5) evolvers: generate a series of atmosphere solutions, of increasing mass
   -`seq1D`: trivial "evolution" since solutions independent
   -`evolve_dM`: take a single step from a static solution to an evolved one
   -`evolve_DM`: take multiple steps over a larger mass range
"""

import sys
import numpy as np
from utils.constants import G, Me, Re, Msun, sigma, Pdisk, Tdisk, RHill, \
                            kdust, AU, Myr, kb, mp
from utils.parameters import FT, FSigma, mstar, rhoc, a
from utils.mynumsci import nextguess, zbrac, mylogspace, ezroot, ngQuad
from scipy.integrate import odeint
from scipy.optimize import brentq, newton, root, fsolve
from scipy.interpolate import interp1d
from collections import namedtuple
from utils.userpath import userpath, aydat
#datfolder = userpath + "/dat/MODELS/RadSGPoly2/" 


mu = 2.35
R = kb / (mp*mu)
p4 = 4*np.pi
mcore  = 5*Me 
rcore = (3*mcore/p4/rhoc)**(1./3)
"""All parameters places in `namedtuple`.  Allows attribute reference like
recarray but no type declatation required! Use namedtuple._replace() to change
values (and create new instance).
"""
Params = namedtuple('Params',
                    'mco, rco, a, R, delad, Po, To,  Mstar, kappa, match')
prms = Params(mcore, rcore, a, R, 2./7, Mstar=mstar*Msun, match='Hill',
              Po=Pdisk(a,mstar, FSigma, FT), To=Tdisk(a,FT), kappa=kdust)

Previ = namedtuple('Previ', 
                   'Mmax, lL, T_m, P_m, r, L, Sold, Mconv, Mcint')
Previunterp = namedtuple('Previunterp', 
                         'Mmax, lL, T, P, m, r, L, Sold, Mconv, Mcint')

def prmsSet(mc, a=a, rhoc=rhoc):
    """create prms with different core mass, radius or disk radius.  
    Other params TODO"""

    rc = (3*mc/p4/rhoc)**(1./3)
    return Params(mc, rc, a, R, 2./7, Mstar=mstar*Msun, match='Hill',
              Po=Pdisk(a,mstar, FSigma, FT), To=Tdisk(a,FT), kappa=kdust)


def evostruct_r(y, r, dt, prev, prms=prms, verbose=0):
    """
    structure eqns. dy/dr for `y` = [m, P , T, L] evolving `dt` from `prev`
    solution

    By setting the global `conv` switch, an inward integration will keep the
    inner regions convective.  Helps keep trial solutions well behaved by
    avoiding spurious radiative zones for too low values of Luminsoity.
    Currently the switch is set when luminosity becomes negative or when
    delrad > 10*delad, whichever happens first.
    
    A numerical concern is that adaptive stepsize integrators could take a big
    trial step and trigger the convective switch too early.  Forcing output at
    regular intervals preventing big steps.  This is done in `error2D`

    Phycially there could be real radiative windows for detailed opacity laws
    that we aren't considering yet.  The L<0 condition avoids this problem, but
    other choices including the delrad switch don't.  It's always possible to
    check if delrad dropped below delad after the fact.

    Not intended for static solutions and will break if dt and prev not given.
    """
    try:
       evostruct_r.conv
    except NameError:
       evostruct_r.conv = 0
    try:
       evostruct_r.Lmin
    except NameError:
       evostruct_r.Lmin = 0
    m, P, T, L = y
    rho = P/prms.R/T
    delad = prms.delad
    if evostruct_r.conv: #convective swith has already been flipped
        Del = delad
    elif L < evostruct_r.Lmin: #Luminosity switch
        evostruct_r.conv = 1 #flip convective switch
        if verbose:
            print "conv. switch flipped at L = ", L
        Del = delad
    else:
        delrad = 3/(16*p4*G*sigma)*prms.kappa(T,P) * P / (m * T**4) * L
        if delrad < delad: #check if convective
            Del = delrad
        else:
            Del = delad
            if delrad > 10*delad:
                evostruct_r.conv = 1 #force to stay convective
                if verbose:
                   print "conv. switch flipped at delrad = ", delrad
    
    dmdr = p4 * rho * r**2
    dPdr = - rho * G * m / r**2
    dTdr = Del * T / P * dPdr
    
    CP = prms.R/delad
    #check if previous solution on known adiabat
    if m < prev.Mconv: 
        Sold = prev.Sold
    else: #TODO: could include S_m in `prev`
        Told = prev.T_m(m)
        Pold = prev.P_m(m)
        Sold = CP*np.log( Told / Pold**delad )
    DS = CP*np.log( T / P**delad ) - Sold
    dLdm = -T * DS / dt
   
    return dmdr, dPdr, dTdr, dLdm*dmdr
    

def polyradstruct_r(y, r, dt = 0, prms = prms, prev = 0, Lfix = 0, vir = 0):
    """
    structure eqns. dy/dr for "y" = [m, P , T, L]  
    """
    m, P, T, L = y
    rho = P/prms.R/T
    delad = prms.delad
    if Lfix != 0:
        L = Lfix
    delrad = 3/(16*p4*G*sigma)*prms.kappa(T,P) * P / (m * T**4) * L 
    Del = min(delrad, delad)
    
    dmdr = p4 * rho * r**2
    dPdr = - rho * G * m / r**2
    dTdr = Del * T / P * dPdr
    if prev == 0:
        dLdm = 0
    else:
        CP = prms.R/delad
        Told = prev.T_m(m)
        Pold = prev.P_m(m)
        DS = CP*np.log(T / Told * (Pold / P)**delad )
        if dt == 0: #energy equation just for tracking purposes.
            dLdm = -T * DS #actually dLdmdot
        else:#actually integrating energy equation
            dLdm = -T * DS / dt
    
    if vir == 0:
        return dmdr, dPdr, dTdr, dLdm*dmdr
    else:#TODO: calc virial terms
        pass
        #return dmdlnP , drdlnP , - Gm / r * dmdlnP , - 3 * R * T * dmdlnP


def polyradstruct_m(y, m, dt = 0, prms = prms):#, prev = 0, Lfix = 0, vir = 0):
    """
    structure eqns. dy/dm for "y" = [r, P , T, L]

    using mass as the dependent variable appears much worse than r, at least for
    shooting inwards.  Presumably the possible divergence of 1/r**2 (continuity)
    and 1/r**4 (HB).  Clearly no such divergence is possible when the range of r
    is controlled as the dependent variable. In principle, should also check
    error behavior of lnr and lnm eqns...
    """
    r, P, T, L = y
    rho = P/prms.R/T
    delad = prms.delad
    
    delrad = 3/(16*p4*G*sigma)*prms.kappa(T,P) * P / (m * T**4) * L 
    Del = min(delrad, delad)
    
    drdm = 1/(p4 * rho * r**2)
    dPdm = - G * m / (p4 * r**4)
    dTdm = Del * T / P * dPdm
    dLdm = 0
    return drdm, dPdm, dTdm, dLdm


def delRad(y, kappa):
    """
    radiative lapse rate

    Parameters
    ----------
    y : 4 element iterable
        `y` = [m, P , T, L]
    kappa : function
        Ross. mean opacity, takes `T` and `P` as arguments
    """
    m , P, T, L = y
    return 3/(16*p4*G*sigma)* kappa(T,P) * P / (m * T**4) * L 


def polyradstruct_r_simp(y, r, dt = 0, prms = prms, prev = 0):
    """
    simple version of structure eqns. dy/dr for "y" = [m, P , T, L].

    Not currently used, but the simplest form of the full time-dependent 
    structure eqns. for polytrope.  Not intended for static solutions and 
    will break if dt and prev not specified.
    """
    m, P, T, L = y
    rho = P/prms.R/T
    delad = prms.delad
    delrad = 3/(16*p4*G*sigma)*prms.kappa(T,P) * P / (m * T**4) * L 
    Del = min(delrad, delad)
    
    dmdr = p4 * rho * r**2
    dPdr = - rho * G * m / r**2
    dTdr = Del * T / P * dPdr
    
    CP = prms.R/delad
    Told = prev.T_m(m)
    Pold = prev.P_m(m)
    DS = CP*np.log(T / Told * (Pold / P)**delad )
    dLdm = -T * DS / dt
   
    return dmdr, dPdr, dTdr, dLdm*dmdr


def rout(Mtot, prms = prms):
    match = prms.match
    if match == "Hill":
        return RHill(Mtot, prms.a , prms.Mstar)
    elif match == "Bondi":
        nfit = 1
    elif match > 0:
        nfit = match
    else:
        sys.exit("bad matching option")
    return nfit * G * Mtot / prms.R / prms.To


def error1D(l10L,Mtot, dt = 0, x = 'r', prms = prms, sargs = (), full_out = 0):
    """
    Computes matching error for core mass.

    Parameters
    ----------
    l10L : float
        Log10 of trial luminosity [ergs/s].  Log must be used to keep guesses
        positive in the bracketing routine.
    Mtot : float
        mass of core and atmosphere out to matching radius
    sargs : tuple , optional
        extra arguments to be passed to structure equations (prev, Lfix, vir)

    Usage:
    ------
    >>> s2.error1D(24.5, 5.8*s2.Me)
    """
    L = 10**l10L
    rmatch = rout(Mtot, prms)
    messg = 0
    if x == 'r':
        f = polyradstruct_r
        bc = [Mtot, prms.Po, prms.To, L]
        xlims = [rmatch, prms.rco]
        bcin = prms.mco
    else:
        f = polyradstruct_m
        bc = [rmatch, prms.Po, prms.To, L]
        xlims = [Mtot, prms.mco]
        bcin = prms.rco
    if not full_out:
        y = odeint(f, bc, xlims, args=(dt, prms)+sargs)
    else:
        y, messg = odeint(f, bc, xlims, args = (dt, prms)+sargs,
                          full_output=True)
    coreint = y[-1,0] #core mass or radius depending
    err = coreint/bcin - 1
    if not full_out:
        return err
    else:
        return err, messg


def errorrange(lL1, lL2, Mtot, x='r', num=50, prms=prms):
    """for visual inspection of the error behavior"""
    lLs = np.linspace(lL1,lL2, num)
    errs = np.zeros(num)
    for i, lL in enumerate(lLs):
        errs[i] = error1D(lL, Mtot, 0, x, prms)
    return errs


def lLtop(Mi, lL1, lL2, dt=0, x='r', prms=prms):
    """
    find luminosity of static solution

    Usage
    -----
    >>> s2.lLtop(5.8*s2.Me,24,25)
    """
    (l10L1, l10L2), success = zbrac(error1D, lL1, lL2, args=(Mi,dt,x,prms))
    if success == False:
        sys.exit('bracket not found')
    return brentq(error1D, l10L1, l10L2, args = (Mi,dt,x,prms))


def minmassfac(prms=prms):
    """
    estimates a minimum total mass (relative to core mass) based on disk mass
    that fills Hill sphere
    """
    rhoo = prms.Po/prms.R/prms.To
    aAU = prms.a*AU
    return 1/(1 - p4*rhoo*aAU**3/9/prms.Mstar)


def seq1D(nM, Mmax, Mmin = 0, follow = 0, lL1 = 25., lL2 = 26., prms=prms,
          savefile = 0):
    """make recarray of solutions.  Works!"""
    Mc = prms.mco
    if Mmin == 0:
        #give 10% safety over uncompressed disk mass might 
        Mmin = minmassfac() * Mc * 1.1 
    Masses = Mc + mylogspace(Mmin - Mc , Mmax - Mc, nM)
    atmser = np.recarray(nM, dtype = [('L', float),('M', float),('R', float),
                                      ('Pc', float),('Tc', float),
                                      ('err', float)])
    if follow:
        y1 , x1 = 0., 0.
    for i, Mi in enumerate(Masses):
        atmser.M[i] = Mi
        atmser.R[i] = rout(Mi, prms)
        if i == 0 or follow == 0:#use fixed guesses
            lLi = lLtop(Mi, lL1 , lL2, prms)
        elif i == 1: #vary about previous
            lLi = lLtop(Mi, 0.9*y1 , 1.1*y1 , prms)
        else:
            yguess = nextguess(Mi, (y0,y1),(x0,x1), maxdy = 10)
            lLi = lLtop(Mi, 0.95*yguess , 1.05*yguess , prms)
        atmser.L[i] = 10**lLi 
        atmser.err[i] = error1D(lLi,Mi, prms = prms)
        #above integration not needed b/c of integration below, but it
        #guarantees that we see the same error that the root finder does.
        y = odeint(polyradstruct_r, [Mi, prms.Po, prms.To, 10**lLi],
                   [atmser.R[i], prms.rco], args=(0, prms))
        atmser.Pc[i] , atmser.Tc[i] = y[-1,1] , y[-1,2]
        if follow: #bookkeeping
            y0, x0 = y1, x1 #previous becomes preprevious
            y1, x1 = lLi, Mi

    if savefile:
        savefile = datfolder + savefile
        #saving as named tuple doesn't seem to work so convert to (ordered)
        #dictionary.
    np.savez_compressed(savefile, params=prms._asdict(), atmser=atmser)
    return atmser


def parammodelsload(file="mc5a10_50.npz"):
    """
    get model parameters and series of models from saved npz file.

    Returns
    -------
    params : namedtuple
    arr : recarray
    """
    npzdat = np.load(datfolder + file)
    #item() takes the dictionary out of the 1 element array
    pdict = npzdat['params'].item() 
    params = Params(**pdict) #turns dictionary back into namedtuple
    arr = npzdat['atmser'].view(np.recarray)

    return params, arr


def profiles(ngrid = 100, file = "mc5a10_50.npz"):
    """return radial profiles
    uses bc's read from file generated by `seq1D`"""
    prms , arr = parammodelsload(file)
    natm = len(arr)
    profs = np.recarray([natm,ngrid],
                        dtype = [('r', float),('m', float),('P', float),
                                 ('TK', float),('L',float),('delrad',float)])
    
    for i, (Mout, Rout, Lout) in enumerate(zip(arr.M , arr.R, arr.L)):
        rgrid = mylogspace(Rout, prms.rco, ngrid)
        y = odeint(polyradstruct_r, [Mout, prms.Po, prms.To, Lout], rgrid,
                   args=(0, prms))
        profs.r[i , :] = rgrid
        profs.m[i , :], profs.P[i , :], profs.TK[i , :], profs.L[i , :] = \
            y[: , 0], y[: , 1], y[: , 2], y[: , 3]

    profs.delrad = delRad([profs.m, profs.P, profs.TK, profs.L], prms.kappa)
    return profs


def sol2prof(prev, prms=prms):
    """turns `makeprev` output into profiles useful for plotting"""
    profs = np.recarray(len(prev.r),
                        dtype = [('r', float), ('m', float), ('P', float),
                                 ('TK', float), ('L',float), ('delrad',float),
                                 ('Mae',float)])
    profs.r = prev.r
    profs.m = prev.T_m.x[::-1]
    profs.P = prev.P_m.y[::-1]
    profs.TK = prev.T_m.y[::-1]
    profs.L = prev.L
    profs.delrad = delRad([profs.m , profs.P , profs.TK , profs.L], prms.kappa)
    profs.Mae = (profs.m - prev.Mcint)/Me
    return profs

def prevunterp(prev):
    """
    remake Previ object without interpolation objects which can't be pickled
    """
    return Previunterp(prev.Mmax, prev.lL, prev.T_m.y, prev.P_m.y, prev.T_m.x
                       , prev.r, prev.L, prev.Sold, prev.Mconv, prev.Mcint)


def prevreterp(prevunt, asdict=1, prms=prms):
    """
    take Previunterp object and re-interpolate T & P
    """
    if asdict:
        prevunt = Previuninterp(**prevunt)
    m, T, P = prevunt.m, prevunt.T, prevunt.P      
    T_m = interp1d(m , T, bounds_error=False, fill_value=prms.To)
    P_m = interp1d(m , P, bounds_error=False, fill_value=prms.Po)
    return Previ(prevunt.Mmax, prevunt.lL, T_m, P_m, prevunt.r, prevunt.L, 
                 prevunt.Sold, prevunt.Mconv, prevunt.Mcint)


def rBondi(m, r, prms):
    """
    find the Bondi radius for a solution. Used for evaluation, not in 
    obtaining solutions.

    When Bondi radius is outside Hill (fit) radius Bondi radius is based on 
    mass within the fit radius.  Bondi radius has little meaning in this case.
    """
    
    if len(m) != len(r):
        sys.exit("mass and radius arrays must be equal length")
    if r[0] > r[-1]: 
        r , m = r[::-1] , m[::-1]
    rBondi_max = G * m[-1] / prms.R / prms.To
    if r[-1] < rBondi_max: #outer radius inside RBondi 
        rBondi = rBondi_max
    else: #compute RBondi
        m_r = interp1d(r , m)
        rBondi = brentq(lambda rr: rr*prms.R*prms.To  - G*m_r(rr),
                        G * prms.mco / prms.R / prms.To, rBondi_max)
    return rBondi


def makeprev(M, lL=0, dt=0, prev=0, prms=prms, ngrid=1000):
    """
    generate data on atmospheric structure in a form that can be used by
    structure eqns. to calculate the T*dS/dt term for luminosity generation.

    Data is output as a named tuple.  Data can be generated for either static
    or evolved solution.  For evolved case, must give both `dt` and `prev`.

    Parameters
    ----------
    M : scalar
        Mass of atmosphere
    lL : scalar , optional for static solutions
        (log10 of) luminosity.
    ngrid : int
        number of grid points for atmosphere solution.  A high number is best
        since 1D interpolation functions.  

    Usage
    -----
    >>> prev = s2.makeprev(M) #static
    >>> prev = s2.makeprev(M, lL, dt, prev)
    """

    Rout = rout(M, prms)
    rgrid = mylogspace(Rout, prms.rco, ngrid)
    if lL == 0:
        lL = lLtop(M, 24, 25, prms=prms)
    if dt == 0:
        y = odeint(polyradstruct_r, [M, prms.Po, prms.To, 10**lL], rgrid,
                   args = (dt, prms))
    else:
        evostruct_r.conv = 0
        y = odeint(evostruct_r , [M, prms.Po, prms.To, 10**lL],
                   rgrid, args = (dt, prev, prms))
    m , P, T , L =  y[: , 0], y[: , 1], y[: , 2], y[: , 3]
    T_m = interp1d(m[::-1] , T[::-1], bounds_error=False, fill_value=prms.To)
    P_m = interp1d(m[::-1] , P[::-1], bounds_error=False, fill_value=prms.Po)
    delad = prms.delad
    CP = prms.R/delad
    Sold = CP*np.log(T[-1] / P[-1]**delad) #use core to define adiabat
    delrad = delRad((m , P, T, L), prms.kappa)
    firstconv = (delrad > delad).argmax()
    #above returns first True, i.e. convective index closest to top.
    #WARNING: if isolated upper convective region exists this will fail (to give
    #         desired result).
    Mconv = m[firstconv]
    return Previ(M , lL, T_m , P_m, rgrid, L, Sold, Mconv, m[-1])


def error2D(x, Mtot, prev, prms=prms, Lnorm=0, full_out=0, dt_log=1,
            num=200):   
    """
    compute matching errors for core mass and luminosity

    Assumed that we want zero Luminosity at core.    
    The global varialble `conv` is set to zero here and then switched within the
    structure equation function, `evostruct_r`

    Parameters
    ----------
    x : two element ndarray, tuple or list
       guess for solution, log10(L) and dt (or log10(dt), the default) [cgs]   
    Mtot : float
        total mass [cgs] of core plus atmosphere to outer boundary
    prev : named tuple 
        previous solution data for structure equations    
    prms : named tuple, optional
        model parameters (core mass, opacity, disk properties, etc.)
    Lnorm : positive float or False
        scaling factor for core luminosity error.  If 0 or 1, then scales to
        outer luminosity.  If < 1e6, then adjusts that standard estimate by
        factor Lnorm (i.e. factors < 1 reduce error).  If > 1e6, then supplied
        value is used as reference luminosity.  This factor helps mean sq.
        errors, (err1**2 + err2**2) have meaningful contributions from both
        terms.
    full_out : int, optional
        If 0, then just returns error.  If 1, then returns ode-int messages.
        If 2 returns atmosphere profile (currently hardcoded to 200 levels).

    Example
    -------
    >>> prev = s2.makeprev(5.8*s2.Me)
    >>> s2.error2D([24.5, .2*s2.Myr], 5.808*s2.Me, prev, dt_log=0)
    """

    #print "core mass {:.2f} Me".format(prms.mco/Me)
    lL = x[0]
    if dt_log:
        dt = 10**x[1]
    else:
        dt = x[1]
    rmatch = rout(Mtot, prms)
    f = evostruct_r
    f.conv, f.Lmin = 0, 1e-4 * 10**lL #flags for structure eqns.
    bc = [Mtot, prms.Po, prms.To, 10**lL]
    rlims = [rmatch, prms.rco]
    rrange = mylogspace(rlims[0], rlims[-1], num)
    stargs = (dt, prev, prms)
    if not full_out:
        y = odeint(f, bc, rrange, args=stargs)
    elif full_out == 1: #give standard odeint messages
        y , messg = odeint(f, bc, rrange, args=stargs, full_output=True)
    else: #return full solution
        y = odeint(f, bc, rrange, args=stargs)
    Mcint , Lcint = y[-1,0] , y[-1,3] #core mass, Luminosity

    #set scale luminosity errors. A numerical Lnorm > 1e6 will pass
    if Lnorm == 0 or Lnorm == 1:
        Lnorm = 10**lL
    elif Lnorm < 1e6: #multiply error by factor Lnorm
        Lnorm = 10**lL/Lnorm
    
    err = (Mcint/prms.mco - 1, Lcint/Lnorm)
    
    #track guesses if lists have been initialized elsewhere
    try:
       #np.append(error2D.guessHist, [x], axis=0)
       error2D.guessHist.append(list(x))
       error2D.errHist.append(list(err))
    except (NameError, AttributeError):
       pass
    
    #return statments as set by output flag
    if not full_out: 
        return err
    elif full_out == 1: #give standard odeint messages
        return err, messg
    else: #return full solution
        return err, rrange, y #y[: , 0], y[: , 1], y[: , 2], y[: , 3]

    
def errormap(Mtot, lL1, lL2, dt1, dt2, prev, prms=prms, nL=10, ndt=10,
             tspace='log'):
    """map error space in 2D"""
    errmap = np.empty((nL,ndt, 2))
    lL = np.linspace(lL1, lL2, nL)
    if tspace == 'lin':
        dt = np.linspace(dt1,dt2,ndt)
    else:
        dt = mylogspace(dt1,dt2,ndt)
    for (i,j) in np.ndindex(nL , ndt):
        errmap[i,j,:] = error2D((lL[i],dt[j]), Mtot, prev, prms, dt_log=0)

    return lL, dt, errmap

def errmapDE(Mtot, lL1, lL2, lDE1, lDE2, prev, prms=prms, nL=10, nDE=10
             ):
    """map error space in 2D, L vs DE = L * dt"""
    errmap = np.empty((nL,nDE, 2))
    lL = np.linspace(lL1, lL2, nL)
    lDE = np.linspace(lDE1, lDE2, nDE)
    for (i,j) in np.ndindex(nL , nDE):
        lDt = lDE[j] - lL[i]
        errmap[i,j,:] = error2D((lL[i], lDt), Mtot, prev, prms, dt_log=1)

    return lL, lDE, errmap


def root2D(x0, Mtot, prev, prms=prms, Lnorm=1, method='ezroot', full_out=0
           , **rkwargs):
    """
    find root x = (lL, log10(dt)) that gives matched atmosphere solution

    Parameters
    ----------
    x0 : two element ndarray, tuple or list
       initial guess for solution (should be close!)
    Mtot : float
        total mass [cgs] of core plus atmosphere to outer boundary
    prms : named tuple, optional
        model parameters (core mass, opacity, disk properties, etc.)
    prev : named tuple
        previous solution data to pass to structure equations
    Lnorm : float, optional
        factor by which to reduce the Luminosity error
    method : string, optional
        root finding method to pass to scipy.optimize.root, currently only
        `hybr` known to work.
    **rkwargs : extra keyword arguments
        passed to root finder

    Example:
        >>> Mold, Mnew = 5.8*s2.Me, 5.808*s2.Me
        >>> prev = s2.makeprev(Mold)
        >>> xg = s2.staticguess(Mnew, prev)
        >>> s2.root2D(xg, Mnew, prev)
    """

    e2args =  (Mtot, prev, prms, Lnorm)
    if method == 'hybr':
        sys.exit("why the heck in hybr")
        return root(error2D, x0, args = e2args, method = method)
    #if method == 'fsolve':
        #equivalent to hybr
        #return fsolve(error2D, x0, args = e2args, diag = (.1,.1*Myr)
        #              , full_output = full_out)
        #return fsolve(error2D, x0, args = e2args, full_output = full_out)
    if method == 'ezroot':
        return ezroot(error2D, x0, args=e2args, **rkwargs)
    


def staticguess(M, prev, dt_guess=0.01*Myr, prms=prms, dt_log=1, statgrid=50,
                full_out=0, check=0):
    """
    make guess for lL and dt (default is log10) based on a static solution

    called by `evolve_dM`.  The main task here is to calculate the `dt` required
    for quasistatic evolution between two (constant L) solutions.  This guess
    for `dt` is used as starting point to search for an evolving solution.

    Usage:
    ------
    >>> prev = s2.makeprev(5.8*s2.Me)
    >>> lLg, dtg = s2.staticguess(5.808*s2.Me, prev, dt_log=0)
    """
    
#obtain static solution with L fixed through atmosphere (set by `args[3]` in
#`odeint` call to `polyradstruct_r` below).
    lL = lLtop(M, prev.lL , .95*prev.lL, prms = prms)
    rgrid = mylogspace(rout(M, prms), prms.rco, statgrid)
    y = odeint(polyradstruct_r, [M, prms.Po, prms.To, 10**lL], rgrid
               , args=(dt_guess, prms, prev, 10**lL))
    Lint =  y[: , 3] #`Lint` returns L + integral of T dS/dt_guess dm/dr

#determine dt value needed to generate assumed (fixed) luminosity. 
    dt = (1 - Lint[-1]/Lint[0]) * dt_guess
    if 0.999 < dt/dt_guess < 1.001:
        sys.exit("numerical errors due to Luminosity runaway. Try finer r grid "
                 "or better dtguess")
    
    if check: #redo integration with calculated dt.  Not unless numerical errors
              #exist b/c with fixed L energy eqn. doesn't affect structure.
        y = odeint(polyradstruct_r, [M, prms.Po, prms.To, 10**lL], rgrid,
                   args=(dt, prms, prev, 10**lL))
        L =  y[:, 3]
        dt = (1 - L[-1]/L[0]) * dt #rescale dt
    else: #this rescaling should be equivalent (barring numerical errors)
        L = Lint[0] + (Lint - Lint[0])/dt

    if (not full_out) and dt_log:
        return lL, np.log10(dt)
    elif not full_out:
        return lL, dt
    else:
        m , P, T  =  y[:, 0], y[:, 1], y[:, 2]
        return (lL, dt), rgrid, m, P, T, L


def evolve_dM(M, dM=0, xg=[], prevroot=0, prms=prms, prevgrid=1000,
              retry=0, full_out=0, Lnorm=1, verbose=1, 
              solver='ezroot', **rkwargs):
    """advance to mass M + dM by K-H contraction
    
    using a specified guess and/or a static solution to make initial
    estimate of L, dt"""

#obtain previous solution
    if not prevroot: #find static solution for initial mass
        prev = makeprev(M, prms=prms, ngrid=prevgrid)
    elif type(prevroot) == float: #log10(L) provided for static solution
        prev = makeprev(M, lL=prevroot, prms=prms, ngrid=prevgrid)
    #use existing prev
    elif 'Previ' in str(type(prevroot)): 
        prev = prevroot
    else: #e.g. supply pre-previous solution for time dependence?
        sys.exit("chosen option for prevroot not supported")

#increment mass
    if dM == 0: #default increase atmosphere mass by 1%
       dM = 0.01*(M - prms.mco)
    M += dM

    eargs = [M, prev, prms, Lnorm] 

#make guess for lL, dt (static guess only if supplied guess missing or bad)
    guess, guessbad = list(xg), 0 #guess = True if xg not empty
    if guess:
        x0 = xg
        err0 = error2D(xg, *eargs)
        guessbad = np.any(np.abs(err0) > 1e-2)
    if (guessbad) or (not guess):
        xg = staticguess(M, prev, prms=prms, statgrid=50)
        xg = (xg[0] + np.log10(1.2), xg[1]) #bump static luminosity guess by 10%
        errg = error2D(xg, *eargs)
        if guessbad: #decide which to keep
            if abs(err0[1]) > abs(errg[1]):
                x0, err0 = xg, errg
                if verbose:
                    print ("choosing static guess for M = {:.3f} M_E"
                           .format(M/Me))
        else:
            x0, err0 = xg, errg             
    if verbose:
        print "{:.2f} Me guess {} has error {}".format(M/Me, x0, err0)
        
#find root and check if improvement
    sol = root2D(x0, *eargs, method=solver, **rkwargs)
    if np.all(np.abs(sol.fun) > np.abs(err0)):
        sol.x, sol.fun = x0, err0

    if max(np.abs(sol.fun)) > 1e-4 and retry:
        if verbose:
            print "trying again, got stuck with error", sol.fun
        Lfac = 1e-3
        eargs[3] = Lfac * Lnorm #rescale luminosity errors
        x0 = (sol.x[0], sol.x[1] + np.log10(1.1) ) #10% bump in dt
        sol2 = root2D(x0, *eargs, method=solver, **rkwargs)
        sol2.fun[1] /= Lfac #need to compare errors on same scale
        if np.any(np.abs(sol2.fun) < np.abs(sol.fun)):
            sol = sol2
            if verbose:
                print "second time a better error", sol.fun

    return sol
    
def evolve_DM(M0, DM, ep=0.01, savedir='test/', savefile=0, restart=0, 
              prms=prms, **dmargs):
    """take a bunch of steps...
    uses ezroot as default solver.

    Returns (suggested)
    -------
    bcarr :
        (M, lL, ldt, err0, err1)
    profs :
        from `sol2prof`
    params :
        for reference and preproducability
    """
    #setup Mass grid and arrays
    Ma0 = M0 - prms.mco
    if 'float' in str(type(DM)): #to catch numpy.float64
        #find steps needed to get precision in atm. mass, rounding up
        nM = 1 + np.ceil(np.log((Ma0+DM)/Ma0) / np.log(1+ep))        
    elif type(DM) == int: #supply num of steps
        nM = DM
        DM = (Ma0 * (1+ep)**(nM-1)) - Ma0
    Ms = prms.mco + mylogspace(Ma0, Ma0+DM, nM)
    dMs = Ms[1:] - Ms[:-1]
   
    bca = np.recarray(nM, dtype = [('M', float),('lL', float),('dt', float),
                                      ('err0', float),('err1', float)])
    bca.M = Ms

    if restart == 0:
        #find initial static solution
        prev = makeprev(Ms[0], prms=prms)
        bca.lL[0], bca.dt[0] = prev.lL, 0
        bca.err0[0], bca.err1[0] = prev.Mcint/prms.mco - 1, 0
        prof0 = sol2prof(prev, prms)
        
        #use the first profile to robustly set dtype and shape
        profss = np.recarray((nM, len(prof0)), dtype=prof0.dtype)
        profss[0, :] = prof0
        skip = 1
    else: #use restart data
        skip = 0
        if type(restart) == str:
            #reconstruct dictionary from file, reterp (TODO)
            pass
        #take first (two) steps using restart info
        solm3, solm2, solm1 = restart['solm3'], restart['solm2'], restart['solm1']
        Mm3, Mm2, Mm1 = restart['Mm3'], restart['Mm2'], restart['Mm1']
        prev, prms = restart['prev'], restart['prms']   

    #loop over remaining steps
    for i, Mnew in enumerate(Ms[skip:], skip):
        #get Mold and dm
        if i == 0:
            dM = Mnew - Mm1
            Mold = Mm1
        else:
            dM = dMs[i - 1]
            Mold = Ms[i - 1]
        
        #guess root from previous solutions
        if restart == 0 and i < 3: #not enough solutions for guess
            xg = []
        elif restart == 0 and i == 3: #2previous
            xg = ngMdot(Mnew, [sol.x, oldersolx], [Ms[2], Ms[1]])
        elif i == 0: #use restart values to make guess
            xg = ngMdot3(Mnew, [solm1, solm2, solm3], [Mm1, Mm2, Mm3])
        elif i == 1: 
            xg = ngMdot3(Mnew, [sol.x, solm1, solm2], [Ms[0], Mm1, Mm2])
        elif i == 2: 
            xg = ngMdot3(Mnew, [sol.x, oldersolx, solm1], [Ms[1], Ms[0], Mm1])
        else: #up and running
            xg = ngMdot3(Mnew, [sol.x, oldersolx, oldestsolx],  Ms[i-3:i][::-1])

        if i > skip:
            if i > skip + 1:
                oldestsolx = oldersolx
            oldersolx = sol.x            
        
        sol = evolve_dM(Mold, dM, xg, prev, prms, **dmargs)

        bca.lL[i], bca.dt[i] = sol.x[0], 10**sol.x[1]
        bca.err0[i], bca.err1[i] = sol.fun[0], sol.fun[1]
        prev = makeprev(Mnew, bca.lL[i], bca.dt[i], prev, prms)
        prof = sol2prof(prev, prms)
        if i == 0: #happens in restart case only
            profss = np.recarray((nM, len(prof)), dtype=prof.dtype)
        profss[i, :] = prof

    #generate restart data
    restart = {'Mm1' : Ms[-1], 'Mm2' : Ms[-2], 'Mm3' : Ms[-3],
               'solm1' : (bca.lL[-1], np.log10(bca.dt[-1])), 
               'solm2' : (bca.lL[-2], np.log10(bca.dt[-2])),
               'solm3' : (bca.lL[-3], np.log10(bca.dt[-3])),
               'prev' : prev, 'prms' : prms}

    #saving
    savedir = aydat + savedir
    if savefile == 1:
        savefile = savedir + 'mc{:.1f}_{:.1f}_{:.1f}_ep{:.2f}'.format(
            prms.mco/Me, M0/Me, (M0+DM)/Me, ep)
    if savefile != 0:
        np.savez_compressed(savefile + '.npz', params=prms._asdict(), 
                            bcs=bca, profs=profss)
        
        restart_sv = restart.copy()
        restart_sv['prev'] = prevunterp(prev)._asdict()
        restart_sv['prms'] = prms._asdict()
        restartfile =  savefile + '_rstrt.npz'
        #output = open(restartfile, 'wb')
        #pickle.dump(restart_sv, output)
        #output.close()
        np.savez_compressed(restartfile, restart=restart_sv)

    return bca, profss, restart

def ngMdot(Mnew, oldsols, oldMs):
    """next guess is not just linear extrapolation since log(dt) is a 
    differential quantity.  Thus if dt was constant for two previous steps,
    if the step size is halved then dt should be halved (not kept constant
    as linear extrapolation would do)."""

    lLguess = nextguess(Mnew, [oldsols[0][0], oldsols[1][0]], oldMs)

    dtold = 10**oldsols[0][1]
    if oldMs[0] < oldMs[1]:
        sys.exit('newest solutions first!')
    Mdot = (oldMs[0] - oldMs[1]) / dtold
    dt = (Mnew - oldMs[0]) / Mdot

    return (lLguess, np.log10(dt))

def ngMdot3(Mnew, oldsols, oldMs):
    """next guess is not just linear extrapolation since log(dt) is a 
    differential quantity.  Thus if dt was constant for two previous steps,
    if the step size is halved then dt should be halved (not kept constant
    as linear extrapolation would do)."""

    if oldMs[0] < oldMs[1]:
        sys.exit('newest solutions first!')
    lLguess = ngQuad(Mnew, [oldsols[0][0], oldsols[1][0], oldsols[2][0]], 
                        oldMs)

    tm3 = 10**oldsols[2][1]
    tm2 = tm3 + 10**oldsols[1][1]
    tm1 = tm2 + 10**oldsols[0][1]
    #print np.log10([tm1, tm2, tm3])

    t = ngQuad(Mnew, [tm1, tm2, tm3], oldMs)
    
    dt = t - tm1

    return (lLguess, np.log10(dt))
