import numpy
from utils.constants import G, kb, mp, Rb, Me, Re, Msun, RH, RHe, sigma, \
     cmperau, RHill, gammafn, mufn, Rfn, Cvfn, kdust, kdust10, Tdisk, Pdisk, \
     paramsEOS, params
from utils.parameters import FT, FSigma, mstar, Y, delad, rhoc
from RadSGRealGas.profiles import atmload as atmloadreal
from RadSGPoly.profiles_poly import atmload as atmloadpoly
from utils.interpolation_functions import interplog10S
from scipy.interpolate import interp1d

def massent_real(Y, a, log10entropy, atmfile):

    prms = paramsEOS(Me, Re, Y, a, Pd = Pdisk(a, mstar, FSigma, FT), \
                 Td = Tdisk(a, FT), kappa = kdust)
    #log10Sref = interplog10S(numpy.log10(Tcb), numpy.log10(Pcb))
    
    atm = atmloadreal(atmfile, prms = prms)
    model = atm[0]
    param = atm[1]
    prof = atm[2]

    Tc = param.Tc
    Pc = param.Pc
    MB = param.MB

    log10S = interplog10S(numpy.log10(Tc), numpy.log10(Pc))

    for i in range(len(log10S) - 1):
        if log10S[i] < log10S[i + 1]:
            break
    log10S = log10S[:i]
    MB = MB[:i]
    
    fm = interp1d(log10S[::-1], MB[::-1])
    Mtot = float(fm(log10entropy))

    return Mtot

def massent_poly(Y, delad, a, mass, atmfile):

    prms = params(Me, Re, a, delad, Y, gamma = gammafn(delad), R = Rfn(Y), \
              Cv = Cvfn(Y, delad), Pd = Pdisk(a, mstar, FSigma, FT), \
              Td = Tdisk(a, FT), kappa = kdust)
    
    atm = atmloadpoly(atmfile, prms = prms)
    model = atm[0]
    param = atm[1]
    prof = atm[2]

    Tcb = param.Tcb
    Pcb = param.Pcb
    MB = param.MB

    fT = interp1d(MB, Tcb)
    fP = interp1d(MB, Pcb)

    Tm = float(fT(mass))
    Pm = float(fP(mass))

    Tref = numpy.array([85.] * len(param))
    Pref = Pcb
    
    log10Sref = interplog10S(numpy.log10(Tref), numpy.log10(Pref))

    dS = (prms.R / prms.delad) * numpy.log((Tm / Tref) * \
                                           (Pref / Pm)**prms.delad)


    log10S = numpy.log10(10**log10Sref + dS) # interplog10S(numpy.log10(Tc), numpy.log10(Pc))

    for i in range(len(log10S) - 1):
        if log10S[i] < log10S[i + 1]:
            break
    log10S = log10S[:i]
    MB = MB[:i]
    
    fS = interp1d(MB, log10S)
    log10S = float(fS(mass))

    return log10S
    
