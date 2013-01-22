


from RadSGPoly.atmseries_poly import mcrit_disktime as mcd

a = 10
exec_file = 0
save = 1

if exec_file == 1:

    Mcrit00 = 0 * numpy.ndarray(shape = (10), dtype = float)
    Mcrit03 = 0 * numpy.ndarray(shape = (10), dtype = float)

    for k in range(10):
        Y = 0.0
        execfile('coolingtime_plots_sg.py')
        Mcrit00[k] = mcd(Mc, t)
        
        Y = 0.3
        execfile('coolingtime_plots_sg.py')
        Mcrit03[k] = mcd(Mc, t)

        a = a + 10
    
au = numpy.linspace(10, 100, 10)
plt.semilogx(au, Mcrit00, '-o', label = r'$\nabla_{\mathrm{ad}}=2/7$, $Y=0.0$')
plt.semilogx(au, Mcrit03, '-rs', label = r'$\nabla_{\mathrm{ad}}=2/7$, $Y=0.3$')
plt.legend()
plt.title(r'Critical core mass vs. distance for disk lifetime of 3 Myrs')
if save == 1:
    plt.savefig(userpath + \
                '/figs/ModelAtmospheres/RadSelfGravPoly/Mcrit_vs_a_3Myrs.pdf')
plt.show()
