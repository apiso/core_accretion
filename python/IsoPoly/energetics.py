
"""check virial equilibrium of existing models"""
from IsoPoly.isopoly import *

from utils.userpath import userpath
from scipy.interpolate import interp1d

datfolder = userpath + "/dat/MODELS/isopoly/"
p4 = 4*np.pi

def modelseriesload(file = "mc5au10.npz"):
    """
    get model parameters and series of models from saved npz file.

    Returns
    -------
    model : recarray
    arr : recarray
    """
    npzdat = np.load(datfolder + file)
    model = npzdat['model'].view(np.recarray)
    arr = npzdat['atmser'].view(np.recarray)

    return model, arr

def virialterms(file = "mc5au10.npz"):    
    """
    Compute terms in virial equation from a model whose boundary conditions have been saved to file

    Returns
    -------
    virialr : recarray
        contains the virial terms Eg, Iu (-3p/rho dm thermal integral) , bottom (neg) and top (pos) boundary terms,
        sum of terms and error in atmospheric matching.
    """
    
    model , arr = modelseriesload(file)
    Tcb , mco , rhoco , R , Po , delad = model.Tcb[0] , model.mco[0] , model.rhoco[0] , model.R[0] , model.Po[0] , model.delad[0]
    ns = len(arr.Pc) 

    virialr = np.recarray(ns, dtype = [('Eg', float),('Iu', float),('bottom', float),('top', float),('Mvir', float),('Rvir', float),('err', float),('vsum',float),('Etot',float)])

    rco = (3 * mco / p4 / rhoco)**(1./3)    
    cargs = (Po , Tcb , mco , rhoco , delad)
    zeta = 3 * delad / (1 - delad)
    
    for i, (Pc,Pcb) in enumerate(zip(arr.Pc, arr.Pcb)):
        virialr.Mvir[i] , virialr.Rvir[i] , virialr.Eg[i] , virialr.Iu[i] , virialr.err[i] = trialint(Pc, Pcb, *cargs, vir = 1)
        virialr.bottom[i] = - p4 * rco**3 * Pc
        virialr.top[i] = p4 * virialr.Rvir[i]**3 * Pcb
        virialr.vsum[i] = -virialr.Eg[i] + virialr.Iu[i] + virialr.bottom[i] + virialr.top[i]
        virialr.Etot[i] = (1 - 1/zeta) * virialr.Eg[i] + (virialr.top[i] + virialr.bottom[i])/zeta

    return virialr

def profiles(ngrid = 100, file = "mc5au10.npz"):    
    """
    Compute profiles of the atmosphere solutions.  

    Returns
    -------
    profiler : recarray
        contains profiles.
    """
    
    model , arr = modelseriesload(file)
    Tcb , mco , rhoco , R , Po , delad = model.Tcb[0] , model.mco[0] , model.rhoco[0] , model.R[0] , model.Po[0] , model.delad[0]
    natm = len(arr.Pc) 

    profiler = np.recarray([natm,ngrid], dtype = [('P', float),('r', float),('m', float),('TK', float)])
    #warning: 'T' is not a good name since arr.T gives the transpose, not the field 'T'
    
    rco = (3 * mco / p4 / rhoco)**(1./3)    
    cargs = (Po , Tcb , mco , rhoco , delad)

    for i , (Pc , Pcb) in enumerate(zip(arr.Pc, arr.Pcb)):
        lnPcb = np.log(Pcb)
        lnPgrid = np.linspace(np.log(Pc) , lnPcb , ngrid)
        profiler.P[i , :] = np.exp(lnPgrid)     
        f , bc = polystructlnP , [mco, rco]
        y = odeint(f , bc , lnPgrid , args = (lnPcb,Tcb,delad,mco))
        profiler.m[i , :] , profiler.r[i , :] = y[: , 0] , y[: , 1]
        profiler.TK[i , :] = Tcb * (profiler.P[i , :] / Pcb)**delad

    return profiler

def explicitLint(T1, T2, m1, m2, DeltaS, order = 'linear'):
    """
    compute Luminosity (/ Mdot) of convective region by explicit integrations of T Ds dm.  Note that Ds is constant between two adiabats.

    TODO: integration range only covers the lower mass solution.  the extra luminosty from the mass that went from radiative to convective zones is probably small but should be calculated for completeness and to check.

    Parameters
    ----------
    T1, T2 : array or interp1d
        temperature profile (values at the m1,m2 arrays or already interpolated over that range)
    m1 , m2 : array
        mass values (monotonically increasing from core)
    Delta S : float
        change in entropy per unit mass
    """
    
    
    DeltaM = m2[-1] - m1[-1] #change in mass between models
                             #print "DM = ", DeltaM/Me

    if DeltaM >= 0:
        mgrid = m1
    else:
        mgrid = m2
        
    #define mass grid
    mmid = (mgrid[:-1] + mgrid[1:])/2 #mass @ cell middle
    dm = mgrid[1:] - mgrid[:-1]

    #interpolate T values
    if "interpolate" in str(type(T2)):
        Tav = (T1(mmid) + T2(mmid))/2
    else:
        fT1_m = interp1d(m1 , T1 , kind = order)
        fT2_m = interp1d(m2 , T2 , kind = order)
        Tav = (fT1_m(mmid) + fT2_m(mmid))/2

    intTdm = np.sum(Tav*dm)
    
    LonMdot = -intTdm*DeltaS/DeltaM
    
    return LonMdot

def explicitLint_cb(Mcbs, Rcbs, Pcbs, Tcb, R = R, delad = 2./7, i = 0):
    """
    estimate the luminosity gain from the material that goes from radiative (isothermal)
    to convective regions (or vice-versa).

    We integrate to the appropriate pressure level in the isothermal atmosphere.  Behavior
    near entropy minimum is funny, and could be sampled at higher reslutions.  (In general a
    resolution check would be a good idea.)
    """
    
    CP = R/delad
    dm = Mcbs[1] - Mcbs[0]
    dmabs = abs(dm)
    if dm >= 0:
        jlow = 0 # index for low mass convective zone
    else:
        jlow = 1
        
    mcb = Mcbs[jlow]
    Pcb = Pcbs[jlow]
    Rcb = Rcbs[jlow]
    Pcb_other = Pcbs[1-jlow]
    y = odeint(isoint_pr_m, [Pcb , Rcb], [mcb , mcb + dmabs/2] , args = (R*Tcb , mcb))
        #dPdm = - G * Mcbs[0] / (p4 * Rcbs[0]**4 )
        #Prad = Pcbs[0] + dPdm * dmabs/2
    Prad = y[1 , 0]
    dS = (-1)**jlow * CP * np.log((Prad / Pcb_other)**delad) #T = Tcb held fixed
                                                             #note sign change for mass decreasing case
        
    rrad = y[1 , 1]
    RB = G*mcb / (R * Tcb)
    Pr_anal = Pcb * np.exp(RB/rrad - RB/Rcb)
    
    #print "{:3} , dm = {: .1e} , dS = {:5.2f} , Pr = {: .1e}, rr/RB = {: .1e} ".format(i , dm/Me, dS/R , Prad, rrad/RB)
    #print "{:3} , dm = {: .1e} , Pr = {: .1e} , Pr_an = {: .1e} , rr/RB = {: .1e} ".format(i , dm/Me, Prad , Pr_anal , rrad/RB)

    LonMdot_cb = Tcb*dS # TODO compute average T (convective region gets a bit hotter)

    return LonMdot_cb

def explicitLarr(pfs, model, order = 'linear'):
    """
    wrapper for explicitLint function to calculate Luminosity (/ Mdot) of convective region by explicit integrations of T Ds dm.  Takes a range of profiles, generates the interpolation fns for T(m) , calculates DS between adiabats and passes the integration to explicitLint.

    Returns
    -------
    LonMdot : ndarray
        Luminosity (/Mdot) values in erg/g
    """

    nL = pfs.shape[0] - 1
    LonMdot = np.empty(nL, dtype = float)
    LonMdotCB = np.empty(nL, dtype = float)
    delad , R = model.delad[0] , model.R[0]
    Tcb = model.Tcb[0]
    
    m1 = pfs.m[0 , :]
    T1 = pfs.TK[0,:]
    fT1 = interp1d(m1 , T1 , kind = order)
    Pcb1 = pfs.P[0 , -1] 
    Rcb1 = pfs.r[0 , -1]
    
    for i in range(nL):
        #print i
        m2, T2 = pfs.m[i+1 , :] , pfs.TK[i+1 , :]
        fT2 = interp1d(m2 , T2 , kind = order)
        Pcb2, Rcb2 = pfs.P[i+1 , -1] , pfs.r[i+1 , -1]
        DS = (R / delad) * np.log( T2[-1] / T1[-1] * (Pcb1 / Pcb2)**delad )
        LonMdot[i] = explicitLint(fT1, fT2, m1, m2 , DS)       
        LonMdotCB[i] = explicitLint_cb((m1[-1],m2[-1]), (Rcb1, Rcb2) , (Pcb1 , Pcb2) , Tcb, R , delad, i)
        #by doing this reassignment the number of recarray calls (and interpolation fn. creations) is cut in half
        m1 , T1 , fT1 = m2 , T2 , fT2
        Pcb1 , Rcb1 = Pcb2 , Rcb2
        
    return LonMdot , LonMdotCB

def isoint_pr_m(y, m, RT, Mcb):

    """
    structure eqns. dy/dm for "y" = [P, r]
    assumes non-self gravitating isothermal layer

    dP/dm = - (G * Mcb)/ (4 * pi * r**4)
    dr/dm = RT /  (4 * pi * * P * r**2)    
    """
    P =y[0]
    r = y[1]
    
    return - (G * Mcb)/ (p4 * r**4) , RT /  (p4 * P * r**2)

def Lglobal(Ecbs, Pcbs, mpfs, rpfs, Tcb , R = R, delad = 2./7, i = 0):
    """
    use global cooling model to get luminosity (over Mdot)

    Parameters
    ----------
    Ecbs , Pcbs : [float, float]
        total energy in convective zone (, pressure at conv. boundary) for two consecutive models
    mpfs, rpfs : [ndarray, ndarray] 
        arrays of mass (and radius) profiles for two consecutive models
    Tcb, R, delad : float
        standard...
    i : int
        index, not used in calculations, for debugging/analysis
    """
    M0 , M1 = mpfs[0][-1],  mpfs[1][-1]
    R0 , R1 = rpfs[0][-1],  rpfs[1][-1]
    
    dM = M1 - M0
    dEdM = (Ecbs[1] - Ecbs[0])/dM
    
    CV = R * (1/delad - 1)
    eacc = - G / 2 * ( M0/R0 + M1/R1 ) + CV * Tcb

    if dM >= 0:
        jlow = 0
    else:
        jlow = 1
        
    mref , rref = [M0,M1][jlow] , [R0,R1][jlow]
    mvals , rvals = mpfs[1-jlow] , rpfs[1-jlow]
    rint = interp1d(mvals , rvals, kind = 'linear')
    rprime = rint(mref)
    dV = (-1)**jlow * p4/3 * (rprime**3 - rref**3)

    work = - 0.5 * (Pcbs[0] + Pcbs[1]) * dV / dM
    
    return -dEdM , eacc , work

def globalLarr(pfs, model, vr, order = 'linear'):
    """
    wrapper for calculating luminosity (over Mdot) from global model for a seris of models.

    Returns
    -------
    dEdM_neg : ndarray
        -dE/dm for convective zone
    eacc : ndarray
        energy per mass of accreted matter
    work : ndarray
        -P Del V/ DM , the surface work contribution
    """

    nL = pfs.shape[0] - 1
    dEdM_neg , eacc , work = np.empty(nL, dtype = float) , np.empty(nL, dtype = float) , np.empty(nL, dtype = float)
    delad , R = model.delad[0] , model.R[0]
    Tcb = model.Tcb[0]
    
    m0 , r0 = pfs.m[0 , :] , pfs.r[0 , :]
    Pcb0 = pfs.P[0 , -1] 
    Ecb0 = vr.Etot[0]
    
    for i in range(nL):
        m1, r1 = pfs.m[i+1 , :] , pfs.r[i+1 , :]
        Pcb1, Ecb1 = pfs.P[i+1 , -1] , vr.Etot[i+1]

        dEdM_neg[i] , eacc[i] , work[i] = Lglobal([Ecb0, Ecb1], [Pcb0,Pcb1], [m0,m1] , [r0,r1] , Tcb, R , delad, i)       
    
        #by doing this reassignment the number of recarray calls is cut in half
        m0 , r0  = m1 , r1 
        Pcb0 , Ecb0 = Pcb1 , Ecb1
        
    return dEdM_neg , eacc , work
