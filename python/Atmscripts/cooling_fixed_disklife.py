

a = numpy.linspace(10, 100, 10)

from RadSGPoly.atmseries_poly import mcrit_disktime as mcd

a = 10

Mcrit = 0 * numpy.ndarray(shape = (10), dtype = float)

for k in range(10):
    execfile('coolingtime_plots_sg.py')
    Mcrit[k] = mcd(Mc, t)
    a = a + 10

    
au = numpy.linspace(10, 100, 10)
plt.semilogx(au, Mcrit, '-o')
plt.show()
