
"""Module to read EOS tables"""

import numpy

def EOSread(Nt, Np, Ncol, filename):
    
    """
    Creates a 3-d array of the EOS table
    Nt = number of T points in the EOS file
    Np = max number of pressure points in the EOS file
    Ncol = number of columns

    returns the 3-d EOS array and the temperature vector
    """

    f=open(filename, 'r')
    EOS=0*numpy.ndarray(shape=(Nt,Np, Ncol), dtype=float)
    T=0*numpy.ndarray(shape=(Nt), dtype=float)
    Npv=0*numpy.ndarray(shape=(Nt), dtype=float)
    for i in range(Nt):
        line=f.readline()
        T[i]=float(line.split()[0])
        npr=float(line.split()[1])
        Npv[i]=npr
        for j in range(int(npr)):
            line=f.readline()
            for k in range(Ncol):
                EOS[i, j, k]=float(line.split()[k])
    f.close()
    return EOS, T, Npv

def EOSnewread(Nt, Np, Ncol, Nppts, filename):

    """Similar to EOSread, was used when creating the extends EOS tables"""
    
    f=open(filename, 'r')
    EOS=0*numpy.ndarray(shape=(Nt,Nppts, Ncol), dtype=float)
    T=0*numpy.ndarray(shape=(Nt), dtype=float)
    Npv=0*numpy.ndarray(shape=(Nt), dtype=float)
    for i in range(Nt):
        line=f.readline()
        T[i]=float(line.split()[0])
        npr=float(line.split()[1])
        Npv[i]=npr
        for j in range(Nppts):
            line=f.readline()
            for k in range(Ncol):
                EOS[i, j, k]=float(line.split()[k])
        for m in range(int(Npv[i])-Nppts):
            f.readline()
    f.close()
    return EOS, T
    

