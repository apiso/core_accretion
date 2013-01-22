
"""

Module to export the profiles of model atmosphere sets into text files,
and to read model atmospheres into arrays from files. Modified version of write_file_new:
uses atmospheres created by bisection method.

Updated version of write_file_new, improved on write_file_new by adding the
automated computation of the minimum entropy.

"""

from utils.parameters import Mc, rc
from utils.table_size import ncol
from utils.disk_model import Td, Pd
import numpy
from atmprofiles_realgas import atmtot

def writeatm(Smin, Smax, npoints, T1, T2, n, filename11, filename12, \
             filename21, filename22):
    """
    Creates a set of model atmosphere parameters and writes them into a text file.

    Smin, Smax = minimum and maxim log(entropy) for the set of atmospheres
    npoints = number of model atmospheres
    T1, T2 = initial temperature guesses

    Output:
        -filename11, filename12: text files with the following columns for each atmosphere
        (11 for the low mass branch, 12 for the high mass branch):
            t = temperature vector (from Tc to Tb)
            r = r(t)
            m = m(t)
            p = p(t)
            rho = rho(t)
        -filename21, filename22: text files with the following set of parameters for each atmosphere:
            (21 for the low mass branch, 22 for the high mass branch):
            log(S)
            M_ad = mass on convective zone in Earth masses
            M_iso = mass of radiative zone in Earth masses
            rb = radius of convective zone in Earth radii
            RB = Bondi radius in Earth radii
            RHill = Hill radius in Earth radii
            Pc = core pressure
            Pcb = pressure at the convective boundary
            Tc = core temperature
            Tcb = temperature at the convective boundary
            Eg = gravitational energy of convective zone
            U = internal energy of convective zone
            Etot = total energy of convective zone
            Egiso = gravitational energy of convective zone
            Uiso = internal energy of convective zone            
            Etotiso = total energy of radiative zone
            vir = check of virial equilibrium
            L = luminosity at the convective boundary

    """
    #Smin = 8.876123 #shoot.mins(n, T1, T2, 10**(-4), 1, Stry, step)[0]
    x = atmtot(Smin, Smax, npoints, n, T1, T2)
    atmprof1=x[0]
    atmset1=x[1]
    atmprof2=x[2]
    atmset2=x[3]
##    atmset1rev = atmset1[::-1]
##    atmset = 0 * numpy.ndarray(shape = (2 * n - 1, ncol), dtype = float)
##    for i in range(n):
##        atmset[i] = atmset1rev[i]
##    for i in range(n, 2 * n - 1):
##        atmset[i] = atmset2[i - n + 1]
    S = numpy.linspace(Smin, Smax, npoints)
    f = open(filename11, 'wb')
    f.write(" Td(K)=%s  " % str(Td))
    f.write(" Pd(dyn cm^-2)=%s  " % str(Pd))
    f.write(" Mc(g)=%s  " % str(Mc))
    f.write(" rc(cm)=%s  " % str(rc))
    f.write("%s\n" % str(npoints))
    f.write("  ")
    for i in range(npoints):
        f.write(" %s  " % str(S[i]))
        f.write("%s\n" % str(numpy.shape(atmprof1)[1]))
        for j in range(numpy.shape(atmprof1)[1]):
            f.write("  ")
            thelist=list(atmprof1[i,j,:])
            for item in thelist:
                f.write("%s  " % item)
            f.write("\n")
        
    f.close()
    f = open(filename12, 'wb')
    f.write(" Td(K)=%s  " % str(Td))
    f.write(" Pd(dyn cm^-2)=%s  " % str(Pd))
    f.write(" Mc(g)=%s  " % str(Mc))
    f.write(" rc(cm)=%s  \n" % str(rc))
    f.write("%s\n" % str(npoints))
    f.write("  ")
    f.write("log(S)       M_conv(Me)      M_rad(Me)       rb(Re)         RB(Re)     \
    RHill(Re)       Pc(dyn cm^-2)        Pcb(dyn cm^-2)         Tc(K)       Tcb(K)      \
    Eg(erg)         U(erg)       Etot(erg)    Egiso(erg)        Uiso(erg)       Etotiso(erg)        \
    vir        L(erg s^-1)    \n" )
    for i in range(numpy.shape(atmset1)[0]):
        f.write("  ")
        thelist = list(atmset1[i,:])
        thelist.insert(0, S[i])
        for item in thelist:
            f.write("%s  " % item)
        f.write("\n")
    f.close()

    f = open(filename21, 'wb')
    f.write(" Td(K)=%s  " % str(Td))
    f.write(" Pd(dyn cm^-2)=%s  " % str(Pd))
    f.write(" Mc(g)=%s  " % str(Mc))
    f.write(" rc(cm)=%s  " % str(rc))
    f.write("%s\n" % str(npoints))
    f.write("  ")
    for i in range(npoints):
        f.write(" %s  " % str(S[i]))
        f.write("%s\n" % str(numpy.shape(atmprof2)[1]))
        for j in range(numpy.shape(atmprof2)[1]):
            f.write("  ")
            thelist=list(atmprof2[i,j,:])
            for item in thelist:
                f.write("%s  " % item)
            f.write("\n")
        
    f.close()
    f = open(filename22, 'wb')
    f.write(" Td(K)=%s  " % str(Td))
    f.write(" Pd(dyn cm^-2)=%s  " % str(Pd))
    f.write(" Mc(g)=%s  " % str(Mc))
    f.write(" rc(cm)=%s  \n" % str(rc))
    f.write("%s\n" % str(npoints))
    f.write("  ")
    f.write("log(S)       M_conv(Me)      M_rad(Me)       rb(Re)         RB(Re)     \
    RHill(Re)       Pc(dyn cm^-2)        Pcb(dyn cm^-2)         Tc(K)       Tcb(K)      \
    Eg(erg)         U(erg)       Etot(erg)    Egiso(erg)        Uiso(erg)       Etotiso(erg)        \
    vir        L(erg s^-1)    \n" )
    for i in range(numpy.shape(atmset2)[0]):
        f.write("  ")
        thelist = list(atmset2[i,:])
        thelist.insert(0, S[i])
        for item in thelist:
            f.write("%s  " % item)
        f.write("\n")
    f.close()




def atmparamread(filename):
    """Reads a set of model atmospheres parameters from a file into an array."""
    f = open(filename, 'r')
    f.readline()
    line = f.readline()
    #Td = float(line.split()[0])
    #Pd = float(line.split()[1])
    #Mc = float(line.split()[2])
    #rc = float(line.split()[3])
    n = int(line.split()[0])
    f.readline()
    atm = 0*numpy.ndarray(shape=(n, ncol), dtype=float)
    S = 0*numpy.ndarray(shape=(n), dtype=float)
    for i in range(n):
        line = f.readline()
        S[i] = float(line.split()[0])
        for j in range(ncol ):
            atm[i, j] = float(line.split()[j+1])
    f.close()
    return atm, S



def atmprofileread(filename):
    """Reads a set of model atmospheres profiles from a file into an array."""
    f = open(filename, 'r')
    line1 = f.readline()
    Nst = int(line1.split()[-1])
    line = f.readline()
    Np = int(line.split()[1])
    atm = 0*numpy.ndarray(shape=(Nst, Np, 5), dtype=float)
    S = 0*numpy.ndarray(shape=(Nst), dtype=float)
    f = open(filename, 'r')
    f.readline()
    for i in range(Nst):
        line = f.readline()
        S[i] = float(line.split()[0])
        for j in range(Np):
            line = f.readline()
            for k in range(numpy.shape(atm)[-1]):
                atm[i, j, k] = float(line.split()[k])
    f.close()
    return atm, S


    

    


