"""
Module to simply examine the behavior self-gravitating polytropes.
The external pressure provided by an idealized isothermal non-self-gravitating layer.
The family of solutions is defined by the binding pressure, P_CB, and matching is required
at R_CB = R_B / [1/nfit + Exp(P_CB / P_disk)]
nfit = 1 is taken for simplicity.

USAGE (assumes userpath/python folder in python path):
-----
In [1]: from IsoPoly import isopoly as ip
In [2]: arr = ip.sequence(nP = 500, mco = 5*Me, savefile = 'mc5au10_N1000')

"""

import numpy as np
from utils.constants import G , Me , Tdisk , Pdisk , Re , Ms , AU
from utils.parameters import R, rhoc , mstar
from scipy.integrate import odeint
from scipy.optimize import fminbound, brentq
from utils.zbrac import zbrac
from utils.myspace import loglinspace
import matplotlib.pyplot as plt
import sys

mc = 5*Me #a default value, the core mass keyword can be set in all relevant functions
a = 10. #not just a default value, while Po & To can be specified in fn calls, disk radius can't
Po = Pdisk(a,mstar,1,1)
To = Tdisk(a, 1)
Mstar = Ms*mstar #just Msun with typical values
#note that gas const (R) i.e. mean mol. wt., must be varied in parameters file (could be changed).

def polystructlnP(y, lnP, lnPcb, Tcb , delad = 2./7 , mco = mc , sg = 1, vir = 0):
    """
    structure eqns. dy/dlnP for "y" = [m, r]
    dm/dlnP = - 4 * pi * r**4 *P / (G * m)
    dr/dlnP = - r**2 * P / (G * m * rho) = - r**2 * R * T / (G * m),
    T = Tcb *  (P / Pcb)**delad

    KEYWORDS:
    sg --- 1 (0) to include (ignore self-grav of atmosphere).  W/o self-grav atm. mass still increments, but is ignored in HB.
    """
    P = np.exp(lnP)
    T = Tcb * (P / np.exp(lnPcb) )**delad
    r = y[1]
    if sg:
        Gm = G * y[0]
    else:
        Gm = G * mco 
    
    dmdlnP = - 4 * np.pi * r**4 *P / Gm
    drdlnP = - r**2 * R * T / Gm
    
    if vir == 0:
        return dmdlnP , drdlnP
    else:
        return dmdlnP , drdlnP , - Gm / r * dmdlnP , - 3 * R * T * dmdlnP


def trialint(Pc, Pcb, Po = Po , Tcb = To, mco = mc, rhoco = rhoc , delad = 2./7 , sg = 1, match = "Bondi", neg = 0, output = 'err', vir = 0):
    """
    Integrate structure equations vs lnP from core to CB and see if CB temp is matched. 
    Returns fractional error in CB temp ratio - 1, which crosses zero at the root.

    Parameters
    ----------
    Pc , Pcb : floats
        trial values for bottom and top pressure
    Po , Tcb : floats
        disk pressure and CB temperature (equal to disk temp for assumed isothermal atm)
    mco , rhoco : floats
        core mass and density
    sg , neg , vir: bool (but int will do)
        whether to: include self-grav in convective zone , return negative error (to find
        maxima with minimization routimes) , calculate and return virial integrals
    output : string or int
        which values to return (overridded if vir == true)
    match : string or int or float
        "Bondi" or "Hill": where to match to disk T and P or else nfit: number of Bondi radii

    EXAMPLE:
    In[34]:  ip.trialint(10**10.3*ip.Po , 10**2*ip.Po , sg = 1)
    Out[34]: 0.047812360976131174
    """
    
    rco = (mco * 3 / (4 * np.pi * rhoco) )**( 1. / 3)
    lnPcb = np.log(Pcb)
    f , bc = polystructlnP , [mco, rco]
    if vir:
        bc.extend([0.,0.])
    
    y = odeint(f, bc, [np.log(Pc) , lnPcb], args = (lnPcb,Tcb,delad,mco,sg,vir)) #integration of the structure equations
    if sg:
        Gm = G * y[-1,0]
    else:
        Gm = G * mco 
    RB = Gm / (R * Tcb)

    if match == "Bondi":
        nfit = 1
    elif match == "Hill":
        RH = (Gm / G / 3 / Mstar)**(1./3) * a * AU
        nfit = RH / RB
    else: #assumes numerical value gives
        nfit = match
  
    rCBpred = RB / (1./nfit + np.log(Pcb/Po))
    
    error = (4 / np.pi) * np.arctan(y[-1 , 1] / rCBpred) - 1 #keeps error term between -1 and 1
    if neg <> 0:
        error = -error
        
    mout , rout = y[-1,0] , y[-1, 1]
    if vir:
        return mout, rout, y[-1, 2] , y[-1, 3] , error
    if (output == 'err' or output ==1):
        return error
    elif output == 2:
        return rout , rCBpred , error 
    else:
        return mout, rout, error
    
def core_vals(Pcb, Tcb = To, mco = mc , rhoco = rhoc , delad = 2./7):
    """computes 

    Returns
    -------
    rco : core radius
    Pcguess : analytic guess of the core pressure.
        Though non-selfgrav. the estimate is fairly robust due to virial equilibrium (and fixed entropy to relate P & T)
    """

    rco = (mco * 3 / (4 * np.pi * rhoco) )**( 1. / 3)
    RBc = G * mco / (R * Tcb) #Bondi radius of the core mass only
    Pcguess = Pcb * (delad * RBc / rco)**(1/delad)
    
    return rco, Pcguess

def trymanyP(lPcb, gmin, gmax, Po = Po , Tcb = To, mco = mc,  sg = 1, rhoco = rhoc , delad = 2./7 , nP = 20, data = 0, output = 1, plotit = 1, match = 1):
    """
    return closeness of fit for a range of core pressure (log spaced)
    optional plotting option

    EXAMPLE:
    >>> ip.trymanyP(2.5,.7,1.5,nP = 200, sg = 1)

    PARAMETERS:
    lPcb --- log10(Pcb/Po), log overpressure of CB relative to disk
    gmin, gmax --- range of core pessure guesses, relative to deep analytic no-self-grav solution

    KEYWORDS:
    output --- (1) error values only (2) r+ , r- , err
    data --- return array of error values
    plotit --- make plot of error values
    """
    
    Pcb = Po * 10**lPcb
    rco, Pcguess = core_vals(Pcb, Tcb , mco  , rhoco  , delad )
    lPs , err = np.linspace( np.log10(gmin*Pcguess) , np.log10(gmax*Pcguess) , nP) , np.zeros(nP)
    if output == 1:
        for i, lPguess in enumerate(lPs):
            err[i] = trialint(10**lPguess , Pcb , Po  , Tcb , mco , rhoco  , delad , sg , match = match)
        if plotit:
            plt.figure(output)
            plt.plot(lPs , err)
        if data:
            return err
    elif output == 2:
        rp , rm = np.zeros(nP) , np.zeros(nP)
        for i, lPguess in enumerate(lPs):
            rp[i] , rm[i] , err[i] = trialint(10**lPguess , Pcb , Po  , Tcb , mco , rhoco  , delad , sg, output = output, match = match)
        if plotit:
            plt.figure(output)
            plt.semilogy(lPs , rp, 'r')
            plt.semilogy(lPs , rm, 'g')
            plt.figure(1)
            plt.plot(lPs , err)
        if data:
            return err
    elif output == 3:
        mout, rout = np.zeros(nP) , np.zeros(nP)
        for i, lPguess in enumerate(lPs):
            mout[i] , rout[i] , err[i] = trialint(10**lPguess , Pcb , Po  , Tcb , mco , rhoco  , delad , sg, output = output, match = match)
        if plotit:
            plt.figure(output)
            plt.semilogy(lPs , mout/Me)
            plt.semilogy(lPs , rout/Re)
            plt.figure(1)
            plt.plot(lPs , err)
        if data:
            return err
        
def checkerr(lPcb, gmin, gmax, Po = Po , Tcb = To, mco = mc, rhoco = rhoc , delad = 2./7 , match = 1, full_out = 1):
    """
    The error in outer radius (R_CB) reported by trialint has a maximum value for self-grav solutions.
    A solution for a given outer Pressure (P_CB) exists only if the max. error > 0.

    Parameters
    ----------
    lPcb --- log10(Pcb/Po), log overpressure of CB relative to disk
    gmin, gmax --- range of core pessure guesses, relative to deep analytic no-self-grav solution
    
    Returns 
    -------
    (if full_out <> 1)
        maxerror --- the maximum error, if negative there is no solution
    (if full_out = 1)
        Popt/Pcguess --- the core pressure (relative to analytic guess) that maximizes the error
        maxerror 
        ierr --- optimization info (0 if converged)
        numfunc --- number of iterations taken

    Example
    -------
    In: ip.checkerr(2.8,.7,1.5)
    Out: (1.1605519239414521, -0.00037058452256022445, 0, 13)
    """
    
    Pcb = Po * 10**lPcb
    rco, Pcguess = core_vals(Pcb, Tcb , mco  , rhoco  , delad )

    ming, maxg = gmin*Pcguess , gmax*Pcguess
    sg , neg = 1 , 1 
    Popt, min_negerr , ierr , numfunc = fminbound(trialint , ming , maxg , args = (Pcb, Po  , Tcb , mco , rhoco , delad , sg ,match, neg ) , full_output = 1 )
    maxerr = -min_negerr

    if full_out:
        return Popt/Pcguess , maxerr , ierr, numfunc
    else:
        return maxerr
    
def minent(lPcb_low = 2, lPcb_high = 3, Po = Po , Tcb = To, mco = mc, rhoco = rhoc , delad = 2./7 , match = 1):
    """
    find minimum entropy (maximum P_CB) solution.
    Note: Will break for very hgh core masses which lack even fully adiabatic atmospheric solutions.
          Since these masses are beyond the transition mass (for 1D atmospheres to be valid) they are
          beyond our scope.

    Returns
    -------
    lPcb_me : float
        max value of log10(Pcb/Po) that gives a matched solution
    """

    #bracket the root
    gmin , gmax , full_out = 0.7 , 1.5 , 0
    argv = (gmin, gmax, Po , Tcb , mco , rhoco , delad , match, full_out)
    lPbrak,b, success,d = zbrac(checkerr, lPcb_low , lPcb_high , factor = 1.2, args = argv)
    if success == False:
        print "min entropy solution not bracketed"
        sys.exit()
    lPcb_me = brentq(checkerr, lPbrak[0], lPbrak[1], args = argv)
    
    return lPcb_me

def makeguess(lPcb, previ, maxguessstep, Po = Po, verbose = 0):
    """
    use previous solutions to make a guess for the next solutions.
    Uses Newton stepping to make a guess for Pc/Pcb, which varies much more smoothly than Pc itself.
    Then tranlates to values for Pc

    Input
    -----
    lPcb : float
        log10(Pcb/Po)
    previ : list
        2 previous solutions and and corresponding x = log10(Pcb/Po) values
    maxguessstep : int or float
        max allowed (dy/dy_previous) 
    """

    Pc0, Pc1 , x0 , x1 = previ
    Pcb0 , Pcb1 , Pcbnew = Po*10**x0 , Po*10**x1 , Po*10**lPcb
    y0 , y1 = Pc0/Pcb0 , Pc1/Pcb1
    dyguess = (y1 - y0)/(x1 - x0) * (lPcb - x1)
    maxdy = maxguessstep*(y1 - y0)
    if abs(dyguess) > abs(maxdy):
        if verbose:
            print "limiting Newton step"
        dyguess = maxdy
    yguess = (y1 + 0.9*dyguess, y1 + 1.1*dyguess)
    Pcguess = (yguess[0]*Pcbnew , yguess[1]*Pcbnew) #no need to involve ndarrays for two elements
    return Pcguess
        
def findpair(lPcb , gmin = 0.7, gmax = 1.5, prev = 0, Po = Po , Tcb = To, mco = mc, rhoco = rhoc , delad = 2./7, match = 1, output = 1, verbose = 0):
    """
    find pair of converged solutions
    TODO: allow to make better use of previous solutions.  A much more efficient way would be to:
       -start with previous two solutions (high and low)
       -from the linear guess for the best step dPc consider steps 0.5 and 2 times dPc (for instance)
       -bracket the root by outward expansion from this guess (need to write zbrak)
       -find exact root

    Parameters
    ----------
    LPcb : float, log_10(Pcb/Po)
    output : int (see Returns options below)

    Returns (incomplete)
    -------
    [Pc1, Pc2] : float list of core pressure solutions    
    """

    maxguessstep = 100 #limits reach of Newton guesses

    Pcb = Po * 10**lPcb
    rco, Pcguess = core_vals(Pcb, Tcb , mco  , rhoco  , delad )
    tiargs = (Pcb, Po , Tcb , mco , rhoco , delad, 1 , match)

    if prev: #linear (Newton) guess for next solution (+/- 10% brackets)
        lowguess = makeguess(lPcb, prev[0], maxguessstep, Po, verbose)
        highguess = makeguess(lPcb, prev[1], maxguessstep, Po, verbose)  
        lowbracs, fs , successl , trials = zbrac(trialint, lowguess[0], lowguess[1], args = tiargs)
        highbracs, fs , successh, trialh = zbrac(trialint, highguess[0], highguess[1], args = tiargs)
        if successl == 0 or successh == 0:
            if verbose:
                print "bracket failure using previous solutions at", lPcb
                print "guesses were", lowguess, highguess, "compared to", Pcguess
            prev = 0 
    if prev == 0:
        #find pressure that lies between high and low Pc solutions, quitting if there are no solutions  
        Popt_norm, maxerr , ierr , numfunc = checkerr(lPcb, gmin, gmax, Po , Tcb , mco , rhoco , delad , match)
        if maxerr < 0:
            print "no solutions"
            sys.exit()
        Popt = Popt_norm*Pcguess
        lowguess = (.99*Popt, Popt)
        highguess = (Popt, 1.01*Popt)
        #bracket two solutions
        lowbracs, fs , successl , trials = zbrac(trialint, lowguess[0], lowguess[1], args = tiargs)
        highbracs, fs , successh, trialh = zbrac(trialint, highguess[0], highguess[1], args = tiargs)
        if successl == 0 or successh == 0:
            print "bracket failure at", lPcb
            sys.exit()

    #find the (bracketed) roots
    lowroot = brentq(trialint, lowbracs[0] , lowbracs[1], args = tiargs)
    highroot = brentq(trialint, highbracs[0] , highbracs[1], args = tiargs)

    if output ==1:
        return [lowroot, highroot]
    else:
        mlow, rlow, errlow = trialint(lowroot, *tiargs, output = 3)
        mhigh, rhigh, errhigh = trialint(highroot, *tiargs, output = 3)
        return (mlow, mhigh) , (rlow,rhigh) , (errlow,errhigh)
    
def sequence(nP = 20, Po = Po , Tcb = To, mco = mc, rhoco = rhoc , delad = 2./7, match = 1, savefile = 0, follow = 0, verbose = 0, nlog = 0):
    """
    generate sequence of atmospheres for a given core mass and all entropy values.  If file provided data saved to
    two (record) arrays, one for basic parameters that apply to all the models and another describing the bc's for the
    series of models.  (Profiles not generated here, see energetics.profiles).  Data easily loaded (with tags) by
    np.load(file).  See plot_ip.py for example of data loading and plotting.

    Inputs (all optional)
    ------
    nP : int
        number of unique surface pressure (entropy) values, output arrays twice this length b/c of high and low mass branches
    match: string or int or float
        where atm matches to disk.  "Hill" for Hill, number for number of Bondi radii
    savefile : 0 or string
        If 0, no data is saved
        To save data to file, specify an identifier, e.g "mc5au10". Path and file extension (.npz) are added automatically. 
    follow : bool, 
        whether to use previous solutions (after first two) to predict initial guesses
    verbose : bool
    nlog: int
        number of log spaced steps to get finer sampling near max Pcb

    Returns
    -------
    atmser : recarray
        solutions as numpy record array of:  core pressure ; mass, radius and pressure @ CB ; and matching error
    file.npz : file (if specified)
        output file of above recarray plus another with basic model params
    """

    ipargs = (Po , Tcb , mco , rhoco , delad, match)
    ipkwargs = {"Po":Po , "Tcb":Tcb , "mco":mco , "rhoco":rhoco , "delad":delad, "match":match}
    safety = 1e-9
    lPcb_max = minent(2.,3.,*ipargs)
    if nlog ==0:
        lPs = np.linspace(0 , lPcb_max/(1 + safety) , nP)
    else:
        ep = loglinspace(safety,1., nP, nlog) #relative distance from max(log(Pcb)), log steps at beginning give
                                              #finer sampling near this maximum.
        lPs = lPcb_max * (1 - ep[::-1]) #recall that zero here means disk pressure
    atmser = np.recarray(2*nP, dtype = [('Pcb', float),('Pc', float),('Mcb', float),('Rcb', float),('err', float)])
    if follow:
        y1l , y1h , x1 = 0 , 0 , 0#initializing the previous solution variables
    for count , lPi in enumerate(lPs[::-1]):#do the enumeration starting from the end
        i = nP - 1 - count #accounting for start from maximum Pcb (min ent)
        if verbose:
            print i
        j = 2*nP - 1 - i #index for high mass branch
        Pcbi = Po*10**lPi
        (atmser.Pcb[i] , atmser.Pcb[j]) = (Pcbi , Pcbi)
        #find Pc
        if count < 2 or follow == 0:
            atmser.Pc[i] , atmser.Pc[j] = findpair(lPi, prev = 0, verbose = verbose, **ipkwargs)
        else:
            prevsol = (( y0l, y1l, x0, x1 ) , ( y0h, y1h, x0, x1 ))
            atmser.Pc[i] , atmser.Pc[j] = findpair(lPi, prev = prevsol, verbose = verbose, **ipkwargs)
        if follow:                #bookkeeping
            y0l , y0h , x0 = y1l , y1h , x1#previous becomes preprevious
            y1l , y1h , x1 = atmser.Pc[i] , atmser.Pc[j] , lPs[i]
        
        #get atm values
        atmser.Mcb[i], atmser.Rcb[i], atmser.err[i] = trialint(atmser.Pc[i], Pcbi, output = 3, **ipkwargs)
        atmser.Mcb[j], atmser.Rcb[j], atmser.err[j] = trialint(atmser.Pc[j], Pcbi, output = 3, **ipkwargs)
        
    if savefile:
        from utils.userpath import userpath
        savefile = userpath+'/dat/MODELS/isopoly/'+savefile
        model = np.array([(Po,Tcb, mco,rhoco,delad,R)] , dtype = [('Po',float),('Tcb',float),('mco',float), ('rhoco',float),('delad',float),('R',float)])
        np.savez_compressed(savefile, model = model, atmser = atmser)

    return atmser
    
def sequenceDM(nP = 500, Po = Po , Tcb = To, mco = mc, rhoco = rhoc , delad = 2./7, match = 1, savefile = 0, follow = 0, verbose = 0):
    """
    (IN PROGRESS)
    take an array and insert intermediate solutions whereever dm/m exceeds some threshold of order 1%
    Why: the large dm/M seems to cause errors in the PdV term, but since this term is small it's not clear its worth the effort.
    How: use insert command `atmser = insert(atmser, i, (Pcb, Pc, etc..))`
    """
    atmser = sequence(np,Po,Tcb,mco,rhoco,delad,match,savefile,follow,verbose)
    
