"""

Modified shooting code that uses log P rather than T as the integration variable. Modified version of
shooting_odeint
"""


from utils.constants import G, kb, mp, Rb, Me, Re, RH, RHe, RHill
from utils.parameters import rc, Mc, Y, beta, R, a
import numpy
import scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate
from utils.table_size import Nt, Np, Ncol
from utils.disk_model import Td, Pd, rhod
import utils.cb_corrections as cb
from utils.interpolation_functions_logp import interpolationt, interpolations, \
     interpolationdelad, interpolationfact, interpolationrho, \
     interpolationu
from utils.interpolation_functions_logt import interpolationp2
from scipy.integrate import odeint
from scipy.interpolate import interp1d


#-----------------------------------------------------------------------------

def Tcore(Sc, T1, T2, n, tol ):

    """
    Implements the shooting method to solve for a model adiabatic atmosphere. It shoots for the core
    temperature for an atmosphere with a given entropy. The atmosphere and 
    disk parameters (Mc, rc, Td, Pd) are taken from the parameters file.

    Input:
        Sc = log10 of desired entropy
        T1 = first try for core temperature
        T2 = second try for core temperature
        n = number of temperature grid points
        tol = allowable tolerance at the boundary

    Output:
        Tc1 = core temperature for the first solution
        Tc2 = core temperature for the second solution

    """
    
    Tb = cb.cbcorr(Sc)[0]  #temperature at top of convective zone
    Pb = cb.cbcorr(Sc)[1]  #pressure at top of convective zone
    deladcb = cb.delad(Tb, Sc) #delad at top of convective zone
    delo = deladcb * ((Tb/Td)**(4-beta) / (Pb/Pd)) #delzero
    delinf = 1. /(4 - beta)

    def g(p):
        """
        Hydrostatic balance in the radiative zone; it reduces to 1/rb - 1/RHill = (R*Td/(G*M))int_Pb^Pd{g(p)} 
        """
        return - (1./p) * (1 + (delo / delinf) * (p / Pd - 1))**(1. /(4 - beta)) 
    pint = integrate.quad(g, Pb, Pd)[0] #ratio = RB/rb
    if (pint > 0):

        step = 0.5 #adjustable step size to decrease the temperature in case the radius
                        #diverges

        Sd=interpolations(numpy.log10(Td), numpy.log10(Pd)) #disk entropy

        if Sc < Sd: #only look for solutions with entropy less than Sd

            def f(x, logp):
                """
                structure eqns. "x" = [ r , m ]
                dr/dlogp = -P/rho * r**2/(G * m),
                dm/dlogp = -4 * pi/(G * m) * p * r**4
                """
                p = numpy.exp(logp)
                rhop = 10**interpolationrho(numpy.log10(p), Sc)
                return numpy.array([ - x[0]**2 * p/ rhop * 1./(G * x[1]), \
                                     - 4 * numpy.pi * p * x[0]**4 / ( G * x[1]) ]) 
                   
                

            def delta(temp):
                """
                Returns the relative matching error for the radius of the convective zone
                """

                pinit = 10**interpolationp2(numpy.log10(temp), Sc) #finds the corresponding core pressure for the
                                                                   #guessed temperature temp
                logp1 = numpy.linspace(numpy.log(pinit), numpy.log(Pb), n) #pressure grid, linear in natural log p                
                y = odeint(f, [rc, Mc], logp1) #integration of the structure equations
                print "t is %s" % temp 
        
                while ( y[:,0][-1] > G * (y[:,1][-1]) / (R * Td) ): #check for divergence of radius
                    Tbig = temp #keep track of Tbig
                    temp = temp - step
                    pinit = 10**interpolationp2(numpy.log10(temp), Sc)
                    logp1 = numpy.linspace(numpy.log(pinit), numpy.log(Pb), n)
                    y = odeint(f, [rc, Mc], logp1)
                    print "t is %s" % temp 
                    #repeat the integration procedure on the new temp; if still diverges, further
                    #decrease temp until radius no longer diverges 
         
                rcb1 = y[:,0][-1] #radius of convective zone
                Mcb1 = y[:,1][-1] #mass of convective zone
                RHill1 = RHill(Mcb1, a) #Hill radius
                rcb = (1. / RHill1 + (R * Td / (G * Mcb1)) * pint)**(-1)
                
                deltab1 = (rcb1 - rcb) / rcb1  #relative error
                return deltab1

            
            """
            Further on we find the two atmosphere solutions for a given entropy; first
            we bracket two intervals that each contain one of the solutions (the search
            is carried out by subdivision in m equal intervals); then, we
            find the right core temperature for each interval through the bisection method.
            """
            
            nb = 2 #number of solutions we are looking for, 2 is assumed below
            xb1 = 0 * numpy.ndarray(shape = (nb), dtype = float) #vector that will contain the
                            #lower bound on core temperatures for each bracket
            xb2 = 0 * numpy.ndarray(shape = (nb), dtype = float) #vector that will contain the
                            #upper bound on core temperatures for each bracket
            delta1 = 0 * numpy.ndarray(shape = (nb), dtype = float) #creating empty arrays
            delta2 = 0 * numpy.ndarray(shape = (nb), dtype = float)
            nbb = 0 #initialize soultion counter
            m = 1.0 #initial number of sub-intervals
            #jmax = numpy.floor(numpy.log2((T2 - T1) / m)) + 1 #number related to the
                    #maximum number of sub-intervals 
            d1 = delta(T1) #error for T1
            d2 = delta(T2) #error for T2    
            if (d1 * d2) > 0: #the two solutions have to be included in the (T1, T2)
                #interval, otherwise we might miss one of the solutions 
                while  ((T2 - T1)/m) > 0.5: 
                    x = T1
                    dx = (T2 - T1) / m
                    deltac = numpy.ndarray(shape = int(m + 1), dtype = float)
                    deltac[0] = delta(x)
                    for i in range(int(m)):
                        x = x + dx
                        deltac[i + 1] = delta(x)
                        print "delta = %s" % deltac[i + 1]
                        if (deltac[i] * deltac[i+1] < 0): #root exists here
                            xb1[nbb] = x - dx
                            xb2[nbb] = x
                            delta1[nbb] = deltac[i]
                            delta2[nbb] = deltac[i + 1]
                            nbb = nbb + 1
                            if (nbb == nb):
                               break
                    # maybe below not needed b/c caught above 
                    if nbb == nb:
                        break
                    m = m * 2
                if (nbb == 0):
                    print "There is no atmosphere solution for this entropy. Choose a higher one."
                    flag = 1
                    return flag
                if (nbb == 1):
                    print "%s is the minimum entropy." % Sc
                else: #start to refine root positions

                    T11 = xb1[0] #(T11, T21) = bracketed interval for solution 1
                    T21 = xb2[0] 
                    T12 = xb1[1] #(T12, T22) = bracketed interval for solution 2
                    T22 = xb2[1]
                    delta11 = delta1[0] #error at T11
                    delta21 = delta2[0] #error at T21
                    delta12 = delta1[1] #error at T12
                    delta22 = delta2[1] #error at T22
                    print "Solution 1 is in the interval %s, %s. Solution 2 is in the interval %s, %s." \
                      % (T11, T21, T12, T22)

                    #Tmid1 = (T11 + T21) / 2 #midpoint of the interval
                    #deltab1 = delta(Tmid1)

                    if (delta11 < 0): #orient the search such that delta > 0 lies at x+dx
                        tbis = T11
                        dx = T21 - T11    
                    else:
                        tbis = T21
                        dx = T11 - T21
                
                    Tmid1 = tbis
                    deltab1 = delta(Tmid1)
                    
                    while (abs(deltab1) > tol):
                        dx = dx / 2
                        Tmid1 = tbis + dx
                        print "T is %s" % Tmid1
                        deltab1 = delta(Tmid1)
                        print "error is %s" % deltab1
                        if (deltab1 < 0): #could multiple by bracket values to avoid sign defs. above
                            tbis = Tmid1


                    #Tmid2 = (T12 + T22) / 2
                    #deltab2 = delta(Tmid2)

                    if (delta21 < 0):
                        tbis = T12
                        dx = T22 - T12
                    else:
                        tbis = T22
                        dx = T12 - T22

                    Tmid2 = tbis
                    deltab2 = delta(Tmid2)

                    while (abs(deltab2) > tol):
                        dx = dx / 2
                        Tmid2 = tbis + dx
                        print "T is %s" % Tmid2
                        deltab2 = delta(Tmid2)
                        print "error is %s" % deltab2
                        if (deltab2 < 0):
                            tbis = Tmid2
                    return Tmid1, Tmid2
            
            else:
                print "Choose different temperature guesses or you will miss a solution!"               

        else:
            print "The entropy is higher than the disk entropy. Choose another value."
    else:
        print "The radius of the convective boundary is larger than the Bondi radius. This is not a \
        physical solution, choose a lower entropy."
    
#-----------------------------------------------------------------------------

def shoot(Sc, T1, T2, n, tol):
    
    """
    Input:
        Sc = log10 of desired entropy
        T1 = first try for core temperature
        T2 = second try for core temperature
        n = number of temperature grid points
        tol = allowable tolerance at the boundary

    Output:
        array of the following arrays:
            t = array of temperatures (from Tc to Tb)
            r = array of r(t), found from hydrostatic balance
            m = array of m(t), found from hydrostatic balance
            p = array of p(t), found from interpolation of EOS tables
            rho = array of rho(t), found from interpolation of EOS tables
            u = array of u(t), found from interpolation of EOS tables
            fact = array of fact(t) = P/(rho * delad), found from interpolation of EOS tables
    """
    
    Tc = Tcore(Sc, T1, T2, n, tol)
    Tc1 = Tc[0]
    Tc2 = Tc[1]
    Pb = cb.cbcorr(Sc)[1]
    
    def ivp(temp):
        pinit = 10**interpolationp2(numpy.log10(temp), Sc)
        logp1 = numpy.linspace(numpy.log(pinit), numpy.log(Pb), n)

        def f(x, logp):
            """
            Returns dy/dlogp where y = (r, M, Eg, Ei, Iu), with Iu the 3p/rho dm integral in the virial theorem
            """
            p = numpy.exp(logp)
            rhop = 10**interpolationrho(numpy.log10(p), Sc) #rho as a function of p for a given S
            up = 10**interpolationu(numpy.log10(p), Sc) #u as a function of p for a given S
            
            return numpy.array([ - x[0]**2 * p/ rhop * 1./(G * x[1]), \
                                - 4 * numpy.pi / ( G * x[1]) * p * x[0]**4, 4 * numpy.pi * p * x[0]**3, \
                                - 4 * numpy.pi / (G * x[1]) * up * p * x[0]**4, \
                                - 12 * numpy.pi / (G * x[1]) * p**2 / rhop * x[0]**4])

        E0 = G * Mc**2 / rc
        y = odeint(f, [rc, Mc, - E0, E0, E0], logp1, mxstep = 50000)

        RHill1 = RHill(y[:, 1][-1], a)
        RB1 = G * (y[:, 1][-1]) / (R * Td)
        rb = y[:,0][-1]
        t = 10**interpolationt(numpy.log10(numpy.exp(logp1)), Sc)
        rho = 10**interpolationrho(numpy.log10((numpy.exp(logp1))), Sc)
        Pc = numpy.exp(logp1)[0]
        rhs = y[:,4][-1] + 4 * numpy.pi * (rc**3 * Pc - rb**3 * Pb)

        return numpy.array([t, y[:, 0], y[:, 1], numpy.exp(logp1), rho, y[:, 3][-1] - E0, \
                            y[:, 2][-1] + E0, (rhs + y[:, 2][-1]) / abs(y[:, 2][-1]), RB1, RHill1])

    return ivp(Tc1), ivp(Tc2)

#--------------------------------------------------------------------------

    


            
    
